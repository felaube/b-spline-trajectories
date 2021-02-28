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

var = 3;
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
    case 2 % Problemático
        x_obs = 2000;
        y_obs = 0;
        z_obs = 900;
        R_obs = 400;    
    case 3
        x_obs = 4500;
        y_obs = 0;
        z_obs = 850;
        R_obs = 200;
    case 4
        x_obs = [3500 4500];
        y_obs = [10 0];
        z_obs = [1000 850];
        R_obs = [200 200];
    case 5
        x_obs = 5000;
        y_obs = 10;
        z_obs = 1000;
        R_obs = 200;
    case 6
        x_obs = 500;
        y_obs = 10;
        z_obs = 1200;
        R_obs = 200;

end

m = length(x_obs);

% Constants

% Gravitational acceleration [m/s²]
g = 9.81;

% Air density (speed is considered to be AES) [kg/m³]
rho = 1.225;

% Safe distance between aircraft and obstacle [m]
safe_distance = 50;

% Aircraft Parameters
var = 0;
switch var
    case 0 % DECODE 2
        % Oswald factor [-]
        e = 0.8;

        % Wing area [m²]
        S = 1.554;

        % Span [m]
        b = 3.9422;

        % Maximum thrust (all operating engines) [N]
        T_max = 130;

        % Aircraft mass [kg]
        mass = 25;

        % Aspect Ratio [-]
        AR = b^2/S;

        % k value [-]
        k = 0.0334;

        % Cd0 value [-]
        C_D_0 = 0.0715;
        
        gamma_max = deg2rad(20);
        gamma_min = deg2rad(-20);
        
        V_min_drag = sqrt(mass*g*2/(rho*S))*((k/C_D_0)^(1/4));
        T_min = 1/2*rho*S*2*C_D_0.*V_min_drag.^2;
end



%% Discretization of the receding horizon trajectory

% OBS: A SOLUÇÃO É EXTREMAMENTE SENSÍVEL AO VALOR DE n
n = 51;
t_h = 200;
t_s = 0.2;

tal = linspace(0, 1, n);
t_array = tal*t_h;

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
R_int(1, 1) =  t_h/7*(t_array(1));
R_int(2, 1) =  -(6*(t_array(1)^7/7 - (5*t_array(1)^6*t_h)/6 + 2*t_array(1)^5*t_h^2 - (5*t_array(1)^4*t_h^3)/2 + (5*t_array(1)^3*t_h^4)/3 - (t_array(1)^2*t_h^5)/2))/t_h^6;
R_int(3, 1) =  (15*(t_array(1)^7/7 - (2*t_array(1)^6*t_h)/3 + (6*t_array(1)^5*t_h^2)/5 - t_array(1)^4*t_h^3 + (t_array(1)^3*t_h^4)/3))/t_h^6;
R_int(4, 1) =  (-20*(t_array(1)^7/7 - (t_array(1)^6*t_h)/2 + (3*t_array(1)^5*t_h^2)/5 - (t_array(1)^4*t_h^3)/4))/t_h^6;
R_int(5, 1) =  (15*t_array(1)^7)/(7*t_h^6) - (5*t_array(1)^6)/t_h^5 + (3*t_array(1)^5)/t_h^4;
R_int(6, 1) =  (-6*(t_array(1)^7/7 - (t_array(1)^6*t_h)/6))/t_h^6;
R_int(7, 1) =  (t_array(1)^7)/(7*t_h^6);

for i = 2:n

    R_int(1, i) = t_h/7*(t_array(i)/t_h - 1)^7 - t_h/7*(t_array(1)/t_h - 1)^7;

    R_int(2, i) = -(6*(t_array(i)^7/7 - (5*t_array(i)^6*t_h)/6 + 2*t_array(i)^5*t_h^2 - (5*t_array(i)^4*t_h^3)/2 + (5*t_array(i)^3*t_h^4)/3 - (t_array(i)^2*t_h^5)/2))/t_h^6;

    R_int(3, i) = (15*(t_array(i)^7/7 - (2*t_array(i)^6*t_h)/3 + (6*t_array(i)^5*t_h^2)/5 - t_array(i)^4*t_h^3 + (t_array(i)^3*t_h^4)/3))/t_h^6;

    R_int(4, i) = (-20*(t_array(i)^7/7 - (t_array(i)^6*t_h)/2 + (3*t_array(i)^5*t_h^2)/5 - (t_array(i)^4*t_h^3)/4))/t_h^6;

    R_int(5, i) = (15*t_array(i)^7)/(7*t_h^6) - (5*t_array(i)^6)/t_h^5 + (3*t_array(i)^5)/t_h^4;

    R_int(6, i) = (-6*(t_array(i)^7/7 - (t_array(i)^6*t_h)/6))/t_h^6;

    R_int(7, i) = (t_array(i)^7)/(7*t_h^6);
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
C_u(4:7) = 0;
%C_u(4:7) = 0;
C_v(4:7) = 0;
C_w(4:7) = 0;


%% Optimization process

ic = [C_u(4:7); C_v(4:7); C_w(4:7)];

t_diff = diff(t_array);

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
                               t_array, t_diff, safe_distance, ...
                               g, rho, e, S, b, mass, AR, k, C_D_0, ...
                               T_max, T_min, gamma_max, gamma_min);

nonlcon = @(x) vehicle_cons(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                            R, R_dot, R_dot_dot, R_int, ...
                            t_diff, x_d, y_d, z_d, psi_d, gamma_d, ...
                            g, rho, e, S, b, T_max, mass, AR, k, C_D_0, T_min);

tic
%optimalWayPoints = fmincon(objective, ic(:), [],[],[],[],[],[], [],opts);
for i = 1 : 1
    optimalWayPoints = fmincon(objective, ic(:), [],[],[],[],[],[],nonlcon,opts);
end
t = toc;

sprintf("Mean execution time was %.8f", t/i)


%% Plot final trajectory

C_u = [C_u(1); C_u(2); C_u(3); optimalWayPoints(1:4)];
C_v = [C_v(1); C_v(2); C_v(3); optimalWayPoints(5:8)];
C_w = [C_w(1); C_w(2); C_w(3); optimalWayPoints(9:12)];

[x, y, z, u, v, w, ...
          u_dot, v_dot, w_dot, ...
          u_dot_dot, v_dot_dot, w_dot_dot, ...
          V, V_dot, ...
          gamma, psi, psi_dot, gamma_dot, ...
          psi_dot_dot, gamma_dot_dot, ...
          phi, n, C_L, C_D, D, T] = plot_trajectory(x_0, y_0, z_0, C_u, C_v, C_w, ...
                                                    R_int, R, R_dot, R_dot_dot, ...
                                                    t_array, t_diff, t_h, ...
                                                    g, rho, S, C_D_0, k, mass, ...
                                                    x_obs, y_obs, z_obs, R_obs, m);


%% Transpose Coefficients Vector just to import more easily into Excel spreadsheet
C_u = C_u';
C_v = C_v';
C_w = C_w';

%% Export workspace variables
save("Workspaces/results")

disp(0);
