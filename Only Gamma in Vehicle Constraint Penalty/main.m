%% Prepare workspace
% close all
clear
% clf

%% User Parameters

% Scaling factors

lambda_p = 1; % position
lambda_s = 1; % speed
lambda_prf = 1; % vehicle contraints
lambda_ob = 1; % obstacles constraints
lambda_t = 1; % terminal attitude
lambda_h = 1; % terminal heading angle
lambda_f = 1; % terminal flight path angle

% Optimization options

%opts = optimoptions('fmincon', 'UseParallel', true);
opts = optimoptions('fmincon');
opts.Display = 'iter';
%opts.Display = 'off';
opts.Algorithm = 'sqp';
opts.MaxFunEvals = inf;
opts.MaxIterations = inf;
%opts.StepTolerance = 1e-16;
%opts.TolX = 1e-16;
%opts.ConstraintTolerance = 1e-03;

% Obstacles

var = 1;

switch var
    case 0
        x_obs = [];
        y_obs = [];
        z_obs = [];
        R_obs = [];
    case 1
        x_obs = 3500;
        y_obs = 10;
        z_obs = 1000;
        R_obs = 200;
    case 2
        x_obs = 2000;
        y_obs = 0;
        z_obs = 900;
        R_obs = 400;
end

m = length(x_obs);


%% Discretization of the receding horizon trajectory

% OBS: A SOLUÇÃO É EXTREMAMENTE SENSÍVEL AO VALOR DE n
n = 40;
t_h = 200;
t_s = 0.2;

tal = linspace(0, 1, n);
t = tal*t_h;

% Bernstein basis functions and derivatives

%{
% Functions in t

R = [(1-t).^6;
     6*t.*(1-t).^5;
     15*t.^2.*(1-t).^4;
     20*t.^3.*(1-t).^3;
     15*t.^4.*(1-t).^2;
     6*t.^5.*(1-t);
     t.^6];

R_dot = [-6*(1-t).^5;
         -6*(6*t-1).*(1-t).^4;
         -30*t.*(3*t-1).*(1-t).^3;
         -60*t.^2.*(2*t-1).*(1-t).^2;
         30*t.^3.*(3*t.^2-5*t + 2);
         6*(5-6*t).*t.^4;
         6*t.^5];

R_dot_dot = [30*(1-t).^4;
             -60*(t-1).^3.*(3*t-1);
             30*(15*t.^2-10*t+1).*(1-t).^2;
             -120*t.*(5*t.^3 - 10*t.^2 + 6*t - 1);
             30*t.^2.*(15*t.^2-20*t+6);
             -60*t.^3.*(3*t-2);
             30*t.^4];
%}

% Functions in tal
%%{
R = [(1-tal).^6;
     6*tal.*(1-tal).^5;
     15*tal.^2.*(1-tal).^4;
     20*tal.^3.*(1-tal).^3;
     15*tal.^4.*(1-tal).^2;
     6*tal.^5.*(1-tal);
     tal.^6];

R_dot = [-6*(1-tal).^5;
         -6*(6*tal-1).*(1-tal).^4;
         -30*tal.*(3*tal-1).*(1-tal).^3;
         -60*tal.^2.*(2*tal-1).*(1-tal).^2;
         30*tal.^3.*(3*tal.^2-5*tal + 2);
         6*(5-6*tal).*tal.^4;
         6*tal.^5];

R_dot_dot = [30*(1-tal).^4;
             -60*(tal-1).^3.*(3*tal-1);
             30*(15*tal.^2-10*tal+1).*(1-tal).^2;
             -120*tal.*(5*tal.^3 - 10*tal.^2 + 6*tal - 1);
             30*tal.^2.*(15*tal.^2-20*tal+6);
             -60*tal.^3.*(3*tal-2);
             30*tal.^4];

%%}

%{
 ##### R_int(1, i) = f(t(i)) - f(t(i-1))

R_int = zeros(7, n);
R_int(1, 1) =  t_h/7*(t(1));
R_int(2, 1) =  -(6*(t(1)^7/7 - (5*t(1)^6*t_h)/6 + 2*t(1)^5*t_h^2 - (5*t(1)^4*t_h^3)/2 + (5*t(1)^3*t_h^4)/3 - (t(1)^2*t_h^5)/2))/t_h^6;
R_int(3, 1) =  (15*(t(1)^7/7 - (2*t(1)^6*t_h)/3 + (6*t(1)^5*t_h^2)/5 - t(1)^4*t_h^3 + (t(1)^3*t_h^4)/3))/t_h^6;
R_int(4, 1) =  (-20*(t(1)^7/7 - (t(1)^6*t_h)/2 + (3*t(1)^5*t_h^2)/5 - (t(1)^4*t_h^3)/4))/t_h^6;
R_int(5, 1) =  (15*t(1)^7)/(7*t_h^6) - (5*t(1)^6)/t_h^5 + (3*t(1)^5)/t_h^4;
R_int(6, 1) =  (-6*(t(1)^7/7 - (t(1)^6*t_h)/6))/t_h^6;
R_int(7, 1) =  (t(1)^7)/(7*t_h^6);


for i = 2:n

    R_int(1, i) = t_h/7*(t(i)) - t_h/7*(t(i-1));

    R_int(2, i) = -(6*(t(i)^7/7 - (5*t(i)^6*t_h)/6 + 2*t(i)^5*t_h^2 - (5*t(i)^4*t_h^3)/2 + (5*t(i)^3*t_h^4)/3 - (t(i)^2*t_h^5)/2))/t_h^6 ...
              +(6*(t(i-1)^7/7 - (5*t(i-1)^6*t_h)/6 + 2*t(i-1)^5*t_h^2 - (5*t(i-1)^4*t_h^3)/2 + (5*t(i-1)^3*t_h^4)/3 - (t(i-1)^2*t_h^5)/2))/t_h^6;

    R_int(3, i) = (15*(t(i)^7/7 - (2*t(i)^6*t_h)/3 + (6*t(i)^5*t_h^2)/5 - t(i)^4*t_h^3 + (t(i)^3*t_h^4)/3))/t_h^6 ...
              - (15*(t(i-1)^7/7 - (2*t(i-1)^6*t_h)/3 + (6*t(i-1)^5*t_h^2)/5 - t(i-1)^4*t_h^3 + (t(i-1)^3*t_h^4)/3))/t_h^6;

    R_int(4, i) = (-20*(t(i)^7/7 - (t(i)^6*t_h)/2 + (3*t(i)^5*t_h^2)/5 - (t(i)^4*t_h^3)/4))/t_h^6 ...
              -(-20*(t(i-1)^7/7 - (t(i-1)^6*t_h)/2 + (3*t(i-1)^5*t_h^2)/5 - (t(i-1)^4*t_h^3)/4))/t_h^6;

    R_int(5, i) = (15*t(i)^7)/(7*t_h^6) - (5*t(i)^6)/t_h^5 + (3*t(i)^5)/t_h^4 ...
              -(15*t(i-1)^7)/(7*t_h^6) - (5*t(i-1)^6)/t_h^5 + (3*t(i-1)^5)/t_h^4;

    R_int(6, i) = (-6*(t(i)^7/7 - (t(i)^6*t_h)/6))/t_h^6 ...
              -(-6*(t(i-1)^7/7 - (t(i-1)^6*t_h)/6))/t_h^6;

    R_int(7, i) = (t(i)^7)/(7*t_h^6) - (t(i-1)^7)/(7*t_h^6);

end
%}

% R_int(1, i) = f(t(i)) - f(t(1))

R_int = zeros(7, n);
R_int(1, 1) =  t_h/7*(t(1));
R_int(2, 1) =  -(6*(t(1)^7/7 - (5*t(1)^6*t_h)/6 + 2*t(1)^5*t_h^2 - (5*t(1)^4*t_h^3)/2 + (5*t(1)^3*t_h^4)/3 - (t(1)^2*t_h^5)/2))/t_h^6;
R_int(3, 1) =  (15*(t(1)^7/7 - (2*t(1)^6*t_h)/3 + (6*t(1)^5*t_h^2)/5 - t(1)^4*t_h^3 + (t(1)^3*t_h^4)/3))/t_h^6;
R_int(4, 1) =  (-20*(t(1)^7/7 - (t(1)^6*t_h)/2 + (3*t(1)^5*t_h^2)/5 - (t(1)^4*t_h^3)/4))/t_h^6;
R_int(5, 1) =  (15*t(1)^7)/(7*t_h^6) - (5*t(1)^6)/t_h^5 + (3*t(1)^5)/t_h^4;
R_int(6, 1) =  (-6*(t(1)^7/7 - (t(1)^6*t_h)/6))/t_h^6;
R_int(7, 1) =  (t(1)^7)/(7*t_h^6);

for i = 2:n

    R_int(1, i) = t_h/7*(t(i)/t_h - 1)^7 - t_h/7*(t(1)/t_h - 1)^7;

    R_int(2, i) = -(6*(t(i)^7/7 - (5*t(i)^6*t_h)/6 + 2*t(i)^5*t_h^2 - (5*t(i)^4*t_h^3)/2 + (5*t(i)^3*t_h^4)/3 - (t(i)^2*t_h^5)/2))/t_h^6;

    R_int(3, i) = (15*(t(i)^7/7 - (2*t(i)^6*t_h)/3 + (6*t(i)^5*t_h^2)/5 - t(i)^4*t_h^3 + (t(i)^3*t_h^4)/3))/t_h^6;

    R_int(4, i) = (-20*(t(i)^7/7 - (t(i)^6*t_h)/2 + (3*t(i)^5*t_h^2)/5 - (t(i)^4*t_h^3)/4))/t_h^6;

    R_int(5, i) = (15*t(i)^7)/(7*t_h^6) - (5*t(i)^6)/t_h^5 + (3*t(i)^5)/t_h^4;

    R_int(6, i) = (-6*(t(i)^7/7 - (t(i)^6*t_h)/6))/t_h^6;

    R_int(7, i) = (t(i)^7)/(7*t_h^6);
end

%R_int(1:7, 1:n) = t_h/7; 

%% Global Trajectory - Demanded Trajectory
% Level Flight
V_d(1:n) = 30;
x_d(1:n) = linspace(0, t_h*V_d(1), n);
y_d(1:n) = 10;
z_d(1:n) = 1000;
psi_d(1:n) = 0;
gamma_d(1:n) = 0;
u_d = V_d.*cos(gamma_d).*cos(psi_d);
v_d = V_d.*cos(gamma_d).*sin(psi_d);
w_d = V_d.*sin(gamma_d);

%% Initial Conditions
u_0 = u_d(1);
v_0 = 0;
w_0 = 0;

u_0_dot = 0;
v_0_dot = 0;
w_0_dot = 0;

u_0_dot_dot = 0;
v_0_dot_dot = 0;
w_0_dot_dot = 0;

x_0 = 0;
y_0 = 10;
z_0 = 1200;

%% Vector of Coefficients
C_u = zeros(7, 1);
C_v = zeros(7, 1);
C_w = zeros(7, 1);

C_u(1) = u_0;
C_v(1) = v_0;
C_w(1) = w_0;

C_u(2) = t_h/6*u_0_dot + C_u(1);
C_v(2) = t_h/6*v_0_dot + C_v(1);
C_w(2) = t_h/6*w_0_dot + C_w(1);

C_u(3) = t_h^2/30*u_0_dot_dot - C_u(1) + 2*C_u(2);
C_v(3) = t_h^2/30*v_0_dot_dot - C_v(1) + 2*C_v(2);
C_w(3) = t_h^2/30*w_0_dot_dot - C_w(1) + 2*C_w(2);

% OBS: A SOLUÇÃO É EXTREMAMENTE SENSÍVEL AO VALOR INICIAL DE C_u(4:7)
% C_u(4:7) = C_u(3);
C_u(4:7) = -10;
%C_u(4:7) = 0;
C_v(4:7) = C_v(3);
C_w(4:7) = C_w(3);


%% Optimization process

ic = [C_u(4:7); C_v(4:7); C_w(4:7)];

t_diff = diff(t);

objective = @(x) cost_function(x, C_u, C_v, C_w, ...
                                 x_0, y_0, z_0, ...
                                 R, R_dot, R_dot_dot, R_int, ...
                                 x_d, y_d, z_d, ...
                                 u_d, v_d, w_d, ...
                                 psi_d, gamma_d, ...
                                 t_h, n, ...
                                 m, x_obs, y_obs, z_obs, R_obs, ...
                                 lambda_p, lambda_s, ...
                                 lambda_prf, lambda_ob, ...
                                 lambda_t, lambda_h, lambda_f, ...
                                 t, t_diff);

nonlcon = @(x) vehicle_cons(x, C_u, C_v, C_w, R, R_dot, R_dot_dot, R_int, t_diff, x_d, y_d, z_d, psi_d, gamma_d);

tic
%optimalWayPoints = fmincon(objective, ic(:), [],[],[],[],[],[], constraint,opts);
optimalWayPoints = fmincon(objective, ic(:), [],[],[],[],[],[],nonlcon,opts);
toc

%% Plot final trajectory

C_u = [C_u(1); C_u(2); C_u(3); optimalWayPoints(1:4)];
C_v = [C_v(1); C_v(2); C_v(3); optimalWayPoints(5:8)];
C_w = [C_w(1); C_w(2); C_w(3); optimalWayPoints(9:12)];

x = x_0 + C_u'*R_int;
y = y_0 + C_v'*R_int;
z = z_0 + C_w'*R_int;
u = C_u'*R;
v = C_v'*R;
w = C_w'*R;
u_dot = C_u'*R_dot/t_h;
v_dot = C_v'*R_dot/t_h;
w_dot = C_w'*R_dot/t_h;
u_dot_dot = C_u'*R_dot_dot/(t_h^2);
v_dot_dot = C_v'*R_dot_dot/(t_h^2);
w_dot_dot = C_w'*R_dot_dot/(t_h^2);

V = sqrt(u.^2 + v.^2 + w.^2);
V_dot = diff(V)./t_diff;

gamma = asin(w./V);
psi = asin(v./(V.*cos(gamma)));

psi_dot = diff(psi)./t_diff;
psi_dot_dot = diff(psi_dot)./t_diff(2:end);

gamma_dot = diff(gamma)./t_diff;
gamma_dot_dot = diff(gamma_dot)./t_diff(2:end);

    %% Constants 
        % Wing area [m²]
    S = 1.554;
    
    % Span [m]
    b = 3.9422;
    
    % Maximum thrust (all operating engines) [N]
    T_max = 130;
    
    % Aircraft mass [kg]
    m = 25;
    
    % Gravitational acceleration [m/s²]
    g = 9.81;
    
    % Air density (speed is considered to be AES) [kg/m³]
    rho = 1.225;
    
    % Aspect Ratio [-]
    %A = b^2/S;
    
    % k value [-]
    k = 0.0334;
    
    % Cd0 value (guessed) [-]
    C_D_0 = 0.0715;
    
    %% Calculate bank angle
    phi = atan2(psi_dot(2:end).*V(3:end).*cos(gamma(3:end)), g*cos(gamma(3:end)) + V(3:end).*gamma_dot_dot);
    
    %% Calculate load factor
    n = (g*cos(gamma(3:end)) + V(3:end).*gamma_dot(2:end))./(g*cos(phi));
    
    %% Calculate thrust
    C_L = 2*n*m*g./(rho*S*V(3:end).^2);
    C_D = C_D_0 + k*C_L.^2;
    
    D = 1/2*rho*S*C_D.*V(3:end).^2;
    
    T = D + m*V_dot(2:end) + m*g*sin(gamma(3:end));
    
    m = length(R_obs);
    
figure(1)
plot(x, y, 'color','k','marker','.','markersize',16)
xlabel("x")
ylabel("y")

figure(2)
plot(x, z, 'color','k','marker','.','markersize',16)
xlabel("x")
ylabel("z")

figure(3)
plot(y, z, 'color','k','marker','.','markersize',16)
xlabel("y")
ylabel("z")

figure(4)
plot(t, u, 'color','k','marker','.','markersize',16)
xlabel("t")
ylabel("u")

figure(5)
plot(t, v, 'color','k','marker','.','markersize',16)
xlabel("t")
ylabel("v")

figure(6)
plot(t, w, 'color','k','marker','.','markersize',16)
xlabel("t")
ylabel("w")

figure(7)
plot(t, V, 'color','k','marker','.','markersize',16)
xlabel("t")
ylabel("V")

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#D9FFFF';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

figure(8)
plot3(x, y, z,'-o','Color','b','MarkerSize',10, 'MarkerFaceColor', color)
plot_obstacle(m, x_obs, y_obs, z_obs, R_obs)
grid on
axis equal
xlabel("x")
ylabel("y")
zlabel("z")

figure(9)
plot(t, z, 'color','k','marker','.','markersize',16)
xlabel("t")
ylabel("z")

figure(10)
plot(t, [1 1 n], 'color','k','marker','.','markersize',16)
xlabel("t")
ylabel("n")

figure(11)
plot(t(3:end), T, 'color','k','marker','.','markersize',16)
xlabel("t")
ylabel("T")


%% Transpose Coefficients Vector just to import more easily into Excel spreadsheet
C_u = C_u';
C_v = C_v';
C_w = C_w';


disp(0);
