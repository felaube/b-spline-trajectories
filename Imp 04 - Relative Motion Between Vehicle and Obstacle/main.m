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
opts.MaxFunEvals = 10000;
%opts.MaxIterations = inf;
%opts.StepTolerance = 1e-16;
%opts.TolX = 1e-16;
%opts.ConstraintTolerance = 1e-03;

% Obstacles

var = 8;
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
    case 1.5
        x_obs = 1000;
        y_obs = 10;
        z_obs = 1000;
        R_obs = 200;
    case 2 % Problemático
        x_obs = 800;
        y_obs = 0;
        z_obs = 900;
        R_obs = 400;    
    case 2.1
        x_obs = [800 1000];
        y_obs = [0 700];
        z_obs = [900 900];
        R_obs = [400 150];   
    case 2.2
        x_obs = [800 1000 3000];
        y_obs = [0 700 500];
        z_obs = [900 900 1000];
        R_obs = [400 150 200];   
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
    case 7
        x_obs = 6000;
        y_obs = 0;
        z_obs = 1000;
        R_obs = 400; 
    case 8
        x_obs = [6000 8000 8750 10000];
        y_obs = [0 -10 10 0];
        z_obs = [1000 1200 900 1000];
        R_obs = [400 200 100 200]; 

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

%% Bernstein basis functions and derivatives

% Functions in tal
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

%% Define current, max and min inputs to avoid local minima:
gamma_current = 0;
psi_current = 0;
V_current = 30;
n_current = 1;
phi_current = 0;

C_L_current = 2*n_current*mass*g./(rho*S*V_current.^2);
C_D_current = C_D_0 + k*C_L_current^2;

D_current = 1/2*rho*S*C_D_current.*V_current.^2;

T_current = D_current + mass*g*sin(gamma_current);

phi_max = deg2rad(40);
phi_min = deg2rad(-40);

n_max = 2;
n_min = -1.5;

%% Global Trajectory - 
% Level Flight
V_g(1:5*n) = 30;
x_g(1:5*n) = linspace(0, 5*t_h*V_g(1), 5*n);
y_g(1:5*n) = 10;
z_g(1:5*n) = 1000;
psi_g(1:5*n) = 0;
gamma_g(1:5*n) = 0;
u_g = V_g.*cos(gamma_g).*cos(psi_g);
v_g = V_g.*cos(gamma_g).*sin(psi_g);
w_g = V_g.*sin(gamma_g);

% Demanded Trajectory
V_d(1:n) = V_g(1:n);
x_d(1:n) = x_g(1:n);
y_d(1:n) = y_g(1:n);
z_d(1:n) = z_g(1:n);
psi_d(1:n) = psi_g(1:n);
gamma_d(1:n) = gamma_g(1:n);
u_d = u_g(1:n);
v_d = v_g(1:n);
w_d = w_g(1:n);

prev_V = [];
prev_x = [];
prev_y = [];
prev_z = [];
prev_psi = [];
prev_gamma = [];
prev_u = [];
prev_v = [];
prev_w = [];
    
%% Initial Conditions
u_0 = u_d(1);
v_0 = v_d(1);
w_0 = w_d(1);

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
    
% OBS: A SOLUÇÃO É EXTREMAMENTE SENSÍVEL AO VALOR INICIAL DE C_u(4:7)
C_u(4:7) = -10;
C_v(4:7) = 0;
C_w(4:7) = 0;
    
for i = 0 : 200
    
    C_u(1) = u_0;
    C_v(1) = v_0;
    C_w(1) = w_0;

    C_u(2) = t_h/6*u_0_dot + C_u(1);
    C_v(2) = t_h/6*v_0_dot + C_v(1);
    C_w(2) = t_h/6*w_0_dot + C_w(1);

    C_u(3) = t_h^2/30*u_0_dot_dot - C_u(1) + 2*C_u(2);
    C_v(3) = t_h^2/30*v_0_dot_dot - C_v(1) + 2*C_v(2);
    C_w(3) = t_h^2/30*w_0_dot_dot - C_w(1) + 2*C_w(2);


    %% Optimization process

    optim_vars = cell(3, 28);
    cost_val(1:28) = Inf;

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
                                   T_max, T_min, gamma_max, gamma_min, ...
                                   3);

    nonlcon = @(x) vehicle_cons(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                R, R_dot, R_dot_dot, R_int, ...
                                t_diff, x_d, y_d, z_d, psi_d, gamma_d, ...
                                g, rho, e, S, b, T_max, mass, AR, k, C_D_0, T_min, ...
                                3);

    tic
    %for i = 1 : 1
        [optimalWayPoints, fval] = fmincon(objective, ic(:), [],[],[],[],[],[],nonlcon,opts);
    %end
    t = toc;

    sprintf("Mean execution time was %.8f", t/i)


    %% Plot final trajectory

    C_u = [C_u(1); C_u(2); C_u(3); optimalWayPoints(1:4)];
    C_v = [C_v(1); C_v(2); C_v(3); optimalWayPoints(5:8)];
    C_w = [C_w(1); C_w(2); C_w(3); optimalWayPoints(9:12)];

    %{
    optim_vars{1, 1} = C_u;
    optim_vars{2, 1} = C_v;
    optim_vars{3, 1} = C_w;
    cost_val(1) = fval;
    %}


    [x, y, z, u, v, w, ...
      u_dot, v_dot, w_dot, ...
      u_dot_dot, v_dot_dot, w_dot_dot, ...
      V, V_dot, ...
      gamma, psi, psi_dot, gamma_dot, ...
      psi_dot_dot, gamma_dot_dot, ...
      phi, load_factor, C_L, C_D, D, T] = plot_trajectory(x_0, y_0, z_0, C_u, C_v, C_w, ...
                                                          R_int, R, R_dot, R_dot_dot, ...
                                                          t_array, t_diff, t_h, ...
                                                          g, rho, S, C_D_0, k, mass, ...
                                                          x_obs, y_obs, z_obs, R_obs, m, ...
                                                          prev_x, prev_y, prev_z);

    u_0 = u(3);
    v_0 = v(3);
    w_0 = w(3);

    u_0_dot = u_dot(3);
    v_0_dot = v_dot(3);
    w_0_dot = w_dot(3);

    u_0_dot_dot = u_dot_dot(3);
    v_0_dot_dot = v_dot_dot(3);
    w_0_dot_dot = w_dot_dot(3);

    x_0 = x(3);
    y_0 = y(3);
    z_0 = z(3);
    
   
    % saving them in previous trajectory
    prev_V = [prev_V V(1:2)];
    prev_x = [prev_x x(1:2)];
    prev_y = [prev_y y(1:2)];
    prev_z = [prev_z z(1:2)];
    prev_psi = [prev_psi psi(1:2)];
    prev_gamma = [prev_gamma gamma(1:2)];
    prev_u = [prev_u u(1:2)];
    prev_v = [prev_v v(1:2)];
    prev_w = [prev_w w(1:2)];
    
    % Global trajectory
    % Remove the first points from the global trajectory
    V_g(1:2) = [];
    x_g(1:2) = [];
    y_g(1:2) = [];
    z_g(1:2) = [];
    psi_g(1:2) = [];
    gamma_g(1:2) = [];
    u_g(1:2) = [];
    v_g(1:2) = [];
    w_g(1:2) = [];

    % Update Demanded Trajectory
    V_d(1:n) = V_g(1:n);
    x_d(1:n) = x_g(1:n);
    y_d(1:n) = y_g(1:n);
    z_d(1:n) = z_g(1:n);
    psi_d(1:n) = psi_g(1:n);
    gamma_d(1:n) = gamma_g(1:n);
    u_d = u_g(1:n);
    v_d = v_g(1:n);
    w_d = w_g(1:n);
    
    pause(0.5)
end

%% Transpose Coefficients Vector just to import more easily into Excel spreadsheet
C_u = C_u';
C_v = C_v';
C_w = C_w';

%% Export workspace variables
save("Workspaces/results")


%{
% Calculate other options to avoid local minima
phi_array = [phi_min, phi_current, phi_max];
T_array = [T_min, T_current, T_max];
n_array = [n_min, n_current, n_max];

calc_states = @(phi, T, n) calculate_states(phi, T, n, ...
                                            gamma_current, psi_current, ...
                                            V_current, D_current, ...
                                            mass, g);                                    

trajectory_index = 2;

for i = 1 : length(phi_array)
    for j = 1 : length(T_array)
        for l = 1 : length(n_array)
            
            % Initial Conditions
            [u_0, v_0, w_0, u_0_dot, v_0_dot, w_0_dot] = calc_states(phi_array(i), T_array(j), n_array(l));

            C_u = zeros(7, 1);
            C_v = zeros(7, 1);
            C_w = zeros(7, 1);

            C_u(1) = u_0;
            C_v(1) = v_0;
            C_w(1) = w_0;

            C_u(2) = t_h/6*u_0_dot + C_u(1);
            C_v(2) = t_h/6*v_0_dot + C_v(1);
            C_w(2) = t_h/6*w_0_dot + C_w(1);

            % OBS: A SOLUÇÃO É EXTREMAMENTE SENSÍVEL AO VALOR INICIAL DE C_u(3:7)
            C_u(3:7) = -10;
            C_v(3:7) = 0;
            C_w(3:7) = 0;
            
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
                                           T_max, T_min, gamma_max, gamma_min, ...
                                           2);

            nonlcon = @(x) vehicle_cons(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                        R, R_dot, R_dot_dot, R_int, ...
                                        t_diff, x_d, y_d, z_d, psi_d, gamma_d, ...
                                        g, rho, e, S, b, T_max, mass, AR, k, C_D_0, T_min, ...
                                        2);    

            ic = [C_u(3:7); C_v(3:7); C_w(3:7)];
            [optimalWayPoints, fval, exitflag, output] = fmincon(objective, ic(:), [],[],[],[],[],[],nonlcon,opts);
            
            %if ((exitflag ~= -2) && (exitflag ~= 2)) && (fval <= 1e10)
            if (exitflag ~= -2) && (fval <= 1e10)
                
                C_u = [C_u(1); C_u(2); optimalWayPoints(1:5)];
                C_v = [C_v(1); C_v(2); optimalWayPoints(6:10)];
                C_w = [C_w(1); C_w(2); optimalWayPoints(11:15)];

               %{ 
                [x, y, z, u, v, w, ...
                u_dot, v_dot, w_dot, ...
                u_dot_dot, v_dot_dot, w_dot_dot, ...
                V, V_dot, ...
                gamma, psi, psi_dot, gamma_dot, ...
                psi_dot_dot, gamma_dot_dot, ...
                phi, load_factor, C_L, C_D, D, T] = plot_trajectory(x_0, y_0, z_0, C_u, C_v, C_w, ...
                                                                    R_int, R, R_dot, R_dot_dot, ...
                                                                    t_array, t_diff, t_h, ...
                                                                    g, rho, S, C_D_0, k, mass, ...
                                                                    x_obs, y_obs, z_obs, R_obs, m);
                %}
                
                optim_vars{1, trajectory_index} = C_u;
                optim_vars{2, trajectory_index} = C_v;
                optim_vars{3, trajectory_index} = C_w;
                cost_val(trajectory_index) = fval;
            end
                                                
            trajectory_index = trajectory_index + 1;
            
        end
    end
end

[~,index] = min(cost_val);

C_u = optim_vars{1, index};
C_v = optim_vars{2, index};
C_w = optim_vars{3, index};

% Define output, objective and constraint functions
% just for debugging purposes
if index == 1
    optimalWayPoints = [C_u(4:7); C_v(4:7); C_w(4:7)];
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
                               T_max, T_min, gamma_max, gamma_min, ...
                               3);

    nonlcon = @(x) vehicle_cons(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                R, R_dot, R_dot_dot, R_int, ...
                                t_diff, x_d, y_d, z_d, psi_d, gamma_d, ...
                                g, rho, e, S, b, T_max, mass, AR, k, C_D_0, T_min, ...
                                3);
else
    optimalWayPoints = [C_u(3:7); C_v(3:7); C_v(3:7)];
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
                                   T_max, T_min, gamma_max, gamma_min, ...
                                   2);

    nonlcon = @(x) vehicle_cons(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                R, R_dot, R_dot_dot, R_int, ...
                                t_diff, x_d, y_d, z_d, psi_d, gamma_d, ...
                                g, rho, e, S, b, T_max, mass, AR, k, C_D_0, T_min, ...
                                2);
end
            
[x, y, z, u, v, w, ...
u_dot, v_dot, w_dot, ...
u_dot_dot, v_dot_dot, w_dot_dot, ...
V, V_dot, ...
gamma, psi, psi_dot, gamma_dot, ...
psi_dot_dot, gamma_dot_dot, ...
phi, load_factor, C_L, C_D, D, T] = plot_trajectory(x_0, y_0, z_0, C_u, C_v, C_w, ...
                                                    R_int, R, R_dot, R_dot_dot, ...
                                                    t_array, t_diff, t_h, ...
                                                    g, rho, S, C_D_0, k, mass, ...
                                                    x_obs, y_obs, z_obs, R_obs, m);
%}