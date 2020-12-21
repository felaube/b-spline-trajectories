%% Prepare workspace
close all
clear

%% User Defined Parameters

% Optimization options

%opts = optimoptions('fmincon', 'UseParallel', true);
opts = optimoptions('fmincon');
%opts.Display = 'iter';
opts.Display = 'off';
opts.Algorithm = 'sqp';
%opts.Algorithm = 'interior-point';
%opts.Algorithm = 'active-set';
%opts.MaxFunEvals = 5000;
%opts.MaxFunEvals = 10000;
%opts.MaxFunEvals = inf;
%opts.MaxIterations = inf;
%opts.StepTolerance = 1e-16;
%opts.TolX = 1e-16;
%opts.ConstraintTolerance = 1e-03;

% Scaling factors
% Default
LAMBDA_P = 1; % position
LAMBDA_S = 10; % speed
LAMBDA_PRF = 10; % vehicle contraints
LAMBDA_OB = 0; % obstacles constraints
LAMBDA_T = 5; % terminal attitude
LAMBDA_H = 1; % terminal heading angle
LAMBDA_F = 1; % terminal flight path angle

% Testing
PARAMETERS = [1, 2, 3, 4;...
              2, 3, 4, 1;...
              3, 4, 1, 2;...
              4, 1, 2, 3;...
              ];
v = 1:10;

PARAMETERS = nchoosek(v,4);
PARAMETERS = [1, 1, 1, 5];
for row = 1 : length(PARAMETERS)
    
    LAMBDA_P = PARAMETERS(row, 1); % position
    LAMBDA_S = PARAMETERS(row, 2); % speed
    LAMBDA_PRF = PARAMETERS(row, 3); % vehicle contraints
    LAMBDA_OB = 0; % obstacles constraints
    LAMBDA_T = PARAMETERS(row, 4); % terminal attitude
    LAMBDA_H = 1; % terminal heading angle
    LAMBDA_F = 1; % terminal flight path angle

    AVOIDING_LOCAL_MINIMA = false;    

    % Discretization of the receding horizon trajectory

    % Number of waypoints in local trajectory
    NUM_LOCAL_WAYPOINTS = 51;
    % Number of traveled points between consecutive steps. 
    NUM_TRAV_WAYPOINTS = 4;
    % Receding horizon time
    HORIZON_TIME = 200;

    tal = linspace(0, 1, NUM_LOCAL_WAYPOINTS);
    t_array = tal*HORIZON_TIME;

    % t_diff, as implemented here, is an array. It could just as well
    % be a constant, since the points are "equitemporal". However,
    % by maintaining it as an array, the script is ready to deal
    % with non "equitemporal" points, if this scenario is ever considered
    % someday
    t_diff = diff(t_array);

    % Obstacles

    % Bear in mind that the method does not behave well when the starting 
    % position of the obstacle is too close to the starting position 
    % of the vehicle
    OBSTACLES_SCENARIO_ID = 8;
    [x_obs, y_obs, z_obs, R_obs, u_obs, v_obs, w_obs, num_obs] = obstacle_scenario(OBSTACLES_SCENARIO_ID);
    
    % Constants

    % Gravitational acceleration [m/s²]
    G = 9.81;

    % Air density (speed is considered to be AES) [kg/m³]
    RHO = 1.225;

    % Minimum separation between aircraft and obstacle [m]
    SAFE_DISTANCE = 100;

    % Aircraft Parameters
    VEHICLE_ID = 1;
    [e, S, b, T_max, mass, AR, k, C_D_0, ...
     gamma_max, gamma_min, phi_max, phi_min, ...
     n_max, n_min, V_min_drag, T_min] = vehicle_parameters(VEHICLE_ID, G, RHO);

    % Global Trajectory
    % Level Flight
    V_g(1:5*NUM_LOCAL_WAYPOINTS) = 30;
    x_g(1:5*NUM_LOCAL_WAYPOINTS) = linspace(0, 5*HORIZON_TIME*V_g(1), 5*NUM_LOCAL_WAYPOINTS);
    y_g(1:5*NUM_LOCAL_WAYPOINTS) = 10;
    z_g(1:5*NUM_LOCAL_WAYPOINTS) = 1000;
    psi_g(1:5*NUM_LOCAL_WAYPOINTS) = 0;
    gamma_g(1:5*NUM_LOCAL_WAYPOINTS) = 0;
    u_g = V_g.*cos(gamma_g).*cos(psi_g);
    v_g = V_g.*cos(gamma_g).*sin(psi_g);
    w_g = V_g.*sin(gamma_g);

    % Demanded Trajectory
    V_d(1:NUM_LOCAL_WAYPOINTS) = V_g(1:NUM_LOCAL_WAYPOINTS);
    x_d(1:NUM_LOCAL_WAYPOINTS) = x_g(1:NUM_LOCAL_WAYPOINTS);
    y_d(1:NUM_LOCAL_WAYPOINTS) = y_g(1:NUM_LOCAL_WAYPOINTS);
    z_d(1:NUM_LOCAL_WAYPOINTS) = z_g(1:NUM_LOCAL_WAYPOINTS);
    psi_d(1:NUM_LOCAL_WAYPOINTS) = psi_g(1:NUM_LOCAL_WAYPOINTS);
    gamma_d(1:NUM_LOCAL_WAYPOINTS) = gamma_g(1:NUM_LOCAL_WAYPOINTS);
    u_d = u_g(1:NUM_LOCAL_WAYPOINTS);
    v_d = v_g(1:NUM_LOCAL_WAYPOINTS);
    w_d = w_g(1:NUM_LOCAL_WAYPOINTS);

    % Initial Conditions
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

    % Vector of Coefficients
    C_u = zeros(7, 1);
    C_v = zeros(7, 1);
    C_w = zeros(7, 1);

    % OBS: A SOLUÇÃO É EXTREMAMENTE SENSÍVEL AO VALOR INICIAL DE C_u(4:7)
    C_u(4:7) = 10;
    C_v(4:7) = 0;
    C_w(4:7) = 0;


    %% Bernstein basis functions and derivatives

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

    R_int = zeros(7, NUM_LOCAL_WAYPOINTS);
    R_int(1, 1) =  HORIZON_TIME/7*(t_array(1));
    R_int(2, 1) =  -(6*(t_array(1)^7/7 - (5*t_array(1)^6*HORIZON_TIME)/6 + 2*t_array(1)^5*HORIZON_TIME^2 - (5*t_array(1)^4*HORIZON_TIME^3)/2 + (5*t_array(1)^3*HORIZON_TIME^4)/3 - (t_array(1)^2*HORIZON_TIME^5)/2))/HORIZON_TIME^6;
    R_int(3, 1) =  (15*(t_array(1)^7/7 - (2*t_array(1)^6*HORIZON_TIME)/3 + (6*t_array(1)^5*HORIZON_TIME^2)/5 - t_array(1)^4*HORIZON_TIME^3 + (t_array(1)^3*HORIZON_TIME^4)/3))/HORIZON_TIME^6;
    R_int(4, 1) =  (-20*(t_array(1)^7/7 - (t_array(1)^6*HORIZON_TIME)/2 + (3*t_array(1)^5*HORIZON_TIME^2)/5 - (t_array(1)^4*HORIZON_TIME^3)/4))/HORIZON_TIME^6;
    R_int(5, 1) =  (15*t_array(1)^7)/(7*HORIZON_TIME^6) - (5*t_array(1)^6)/HORIZON_TIME^5 + (3*t_array(1)^5)/HORIZON_TIME^4;
    R_int(6, 1) =  (-6*(t_array(1)^7/7 - (t_array(1)^6*HORIZON_TIME)/6))/HORIZON_TIME^6;
    R_int(7, 1) =  (t_array(1)^7)/(7*HORIZON_TIME^6);

    for i = 2:NUM_LOCAL_WAYPOINTS
        R_int(1, i) = HORIZON_TIME/7*(t_array(i)/HORIZON_TIME - 1)^7 - HORIZON_TIME/7*(t_array(1)/HORIZON_TIME - 1)^7;

        R_int(2, i) = -(6*(t_array(i)^7/7 - (5*t_array(i)^6*HORIZON_TIME)/6 + 2*t_array(i)^5*HORIZON_TIME^2 - (5*t_array(i)^4*HORIZON_TIME^3)/2 + (5*t_array(i)^3*HORIZON_TIME^4)/3 - (t_array(i)^2*HORIZON_TIME^5)/2))/HORIZON_TIME^6;

        R_int(3, i) = (15*(t_array(i)^7/7 - (2*t_array(i)^6*HORIZON_TIME)/3 + (6*t_array(i)^5*HORIZON_TIME^2)/5 - t_array(i)^4*HORIZON_TIME^3 + (t_array(i)^3*HORIZON_TIME^4)/3))/HORIZON_TIME^6;

        R_int(4, i) = (-20*(t_array(i)^7/7 - (t_array(i)^6*HORIZON_TIME)/2 + (3*t_array(i)^5*HORIZON_TIME^2)/5 - (t_array(i)^4*HORIZON_TIME^3)/4))/HORIZON_TIME^6;

        R_int(5, i) = (15*t_array(i)^7)/(7*HORIZON_TIME^6) - (5*t_array(i)^6)/HORIZON_TIME^5 + (3*t_array(i)^5)/HORIZON_TIME^4;

        R_int(6, i) = (-6*(t_array(i)^7/7 - (t_array(i)^6*HORIZON_TIME)/6))/HORIZON_TIME^6;

        R_int(7, i) = (t_array(i)^7)/(7*HORIZON_TIME^6);
    end

    % Define initial states
    gamma_current = gamma_d(1);
    psi_current = psi_d(1);
    V_current = sqrt(u_0^2 + v_0^2 + w_0^2);
    n_current = 1;
    phi_current = 0;

    C_L_current = 2*n_current*mass*G./(RHO*S*V_current.^2);
    C_D_current = C_D_0 + k*C_L_current^2;

    D_current = 1/2*RHO*S*C_D_current.*V_current.^2;

    T_current = D_current + mass*G*sin(gamma_current);

    % Initialize arrays that will store the historical path and states 
    past_V = [];
    past_x = [];
    past_y = [];
    past_z = [];
    past_psi = [];
    past_gamma = [];
    past_u = [];
    past_v = [];
    past_w = [];

    % Initialize variables used to store the objective function value
    % and the optmization results
    optim_vars = cell(3, 28);
    cost_val(1:28) = Inf;

    while true
        % tic
        cost_val(:) = Inf;

        C_u(1) = u_0;
        C_v(1) = v_0;
        C_w(1) = w_0;

        C_u(2) = HORIZON_TIME/6*u_0_dot + C_u(1);
        C_v(2) = HORIZON_TIME/6*v_0_dot + C_v(1);
        C_w(2) = HORIZON_TIME/6*w_0_dot + C_w(1);

        C_u(3) = HORIZON_TIME^2/30*u_0_dot_dot - C_u(1) + 2*C_u(2);
        C_v(3) = HORIZON_TIME^2/30*v_0_dot_dot - C_v(1) + 2*C_v(2);
        C_w(3) = HORIZON_TIME^2/30*w_0_dot_dot - C_w(1) + 2*C_w(2);


        %% Optimization process

        ic = [C_u(4:7); C_v(4:7); C_w(4:7)];

        objective = @(x) cost_function(x, C_u, C_v, C_w, ...
                                       x_0, y_0, z_0, ...
                                       R, R_dot, R_dot_dot, R_int, ...
                                       x_d, y_d, z_d, ...
                                       u_d, v_d, w_d, ...
                                       psi_d, gamma_d, ...
                                       HORIZON_TIME, NUM_LOCAL_WAYPOINTS, ...
                                       num_obs, x_obs, y_obs, z_obs, R_obs, ...
                                       u_obs, v_obs, w_obs, ...
                                       LAMBDA_P, LAMBDA_S, ...
                                       LAMBDA_PRF, LAMBDA_OB, ...
                                       LAMBDA_T, LAMBDA_H, LAMBDA_F, ...
                                       t_array, t_diff, SAFE_DISTANCE, ...
                                       G, RHO, e, S, b, mass, AR, k, C_D_0, ...
                                       T_max, T_min, gamma_max, gamma_min, ...
                                       phi_max, phi_min, ...
                                       3);

        nonlcon = @(x) problem_constraints(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                           x_obs, y_obs, z_obs, R_obs, ...
                                           u_obs, v_obs, w_obs, ...
                                           R, R_dot, R_dot_dot, R_int, ...
                                           t_array, t_diff, NUM_LOCAL_WAYPOINTS, num_obs, ...
                                           x_d, y_d, z_d, psi_d, gamma_d, ...
                                           G, RHO, e, S, b, T_max, mass, AR, k, C_D_0, T_min, ...
                                           SAFE_DISTANCE, 3);

        %tic
        [optimalWayPoints, fval, ~, output] = fmincon(objective, ic(:), [],[],[],[],[],[],nonlcon,opts);
        %t = toc;

        % sprintf("Iterations: %i funcCount: %i", output.iterations, output.funcCount)
        % sprintf("Execution time: %.8f", t)

        optim_vars{1, 1} = [C_u(1); C_u(2); C_u(3); optimalWayPoints(1:4)];
        optim_vars{2, 1} = [C_v(1); C_v(2); C_v(3); optimalWayPoints(5:8)];
        optim_vars{3, 1} = [C_w(1); C_w(2); C_w(3); optimalWayPoints(9:12)];
        cost_val(1) = fval;

        %% Avoiding Minima
        if AVOIDING_LOCAL_MINIMA
            % Calculate other options to avoid local minima
            phi_array = [phi_min, phi_current, phi_max];
            T_array = [T_min, T_current, T_max];
            n_array = [n_min, n_current, n_max];

            calc_states = @(phi, T, NUM_LOCAL_WAYPOINTS) calculate_states(phi, T, NUM_LOCAL_WAYPOINTS, ...
                                                        gamma_current, psi_current, ...
                                                        V_current, D_current, ...
                                                        mass, G);                                    

            trajectory_index = 2;

            % tic
            for i = 1 : length(phi_array)
                for j = 1 : length(T_array)
                    for l = 1 : length(n_array)

                        % Calculate Initial Conditions
                        [u_0, v_0, w_0, u_0_dot, v_0_dot, w_0_dot] = calc_states(phi_array(i), T_array(j), n_array(l));

                        C_u(1) = u_0;
                        C_v(1) = v_0;
                        C_w(1) = w_0;

                        C_u(2) = HORIZON_TIME/6*u_0_dot + C_u(1);
                        C_v(2) = HORIZON_TIME/6*v_0_dot + C_v(1);
                        C_w(2) = HORIZON_TIME/6*w_0_dot + C_w(1);

                        objective = @(x) cost_function(x, C_u, C_v, C_w, ...
                                                       x_0, y_0, z_0, ...
                                                       R, R_dot, R_dot_dot, R_int, ...
                                                       x_d, y_d, z_d, ...
                                                       u_d, v_d, w_d, ...
                                                       psi_d, gamma_d, ...
                                                       HORIZON_TIME, NUM_LOCAL_WAYPOINTS, ...
                                                       num_obs, x_obs, y_obs, z_obs, R_obs, ...
                                                       u_obs, v_obs, w_obs, ...
                                                       LAMBDA_P, LAMBDA_S, ...
                                                       LAMBDA_PRF, LAMBDA_OB, ...
                                                       LAMBDA_T, LAMBDA_H, LAMBDA_F, ...
                                                       t_array, t_diff, SAFE_DISTANCE, ...
                                                       G, RHO, e, S, b, mass, AR, k, C_D_0, ...
                                                       T_max, T_min, gamma_max, gamma_min, ...
                                                       phi_max, phi_min, ...
                                                       2);

                        nonlcon = @(x) problem_constraints(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                                           x_obs, y_obs, z_obs, R_obs, ...
                                                           u_obs, v_obs, w_obs, ...
                                                           R, R_dot, R_dot_dot, R_int, ...
                                                           t_array, t_diff, NUM_LOCAL_WAYPOINTS, num_obs, ...
                                                           x_d, y_d, z_d, psi_d, gamma_d, ...
                                                           G, RHO, e, S, b, T_max, mass, AR, k, C_D_0, T_min, ...
                                                           SAFE_DISTANCE, 2);    

                        ic = [C_u(3:7); C_v(3:7); C_w(3:7)];
                        [optimalWayPoints, fval, exitflag, output] = fmincon(objective, ic(:), [],[],[],[],[],[],nonlcon,opts);

                        if (exitflag >= 0) && (fval <= 1e10)

                            C_u = [C_u(1); C_u(2); optimalWayPoints(1:5)];
                            C_v = [C_v(1); C_v(2); optimalWayPoints(6:10)];
                            C_w = [C_w(1); C_w(2); optimalWayPoints(11:15)];

                            optim_vars{1, trajectory_index} = C_u;
                            optim_vars{2, trajectory_index} = C_v;
                            optim_vars{3, trajectory_index} = C_w;
                            cost_val(trajectory_index) = fval;

                            sprintf("Iterations: %i funcCount: %i", output.iterations, output.funcCount)
                        end

                        trajectory_index = trajectory_index + 1;

                    end
                end
            end
        end
        % t = toc;

        % sprintf("Avoiding minima execution time: %.8f", t)    

        [~,index] = min(cost_val);

        C_u = optim_vars{1, index};
        C_v = optim_vars{2, index};
        C_w = optim_vars{3, index};

        % Redefine objective and constraint functions and fmincon output
        % just for debugging purposes
        if index == 1
            optimalWayPoints = [C_u(4:7); C_v(4:7); C_w(4:7)];
            objective = @(x) cost_function(x, C_u, C_v, C_w, ...
                                           x_0, y_0, z_0, ...
                                           R, R_dot, R_dot_dot, R_int, ...
                                           x_d, y_d, z_d, ...
                                           u_d, v_d, w_d, ...
                                           psi_d, gamma_d, ...
                                           HORIZON_TIME, NUM_LOCAL_WAYPOINTS, ...
                                           num_obs, x_obs, y_obs, z_obs, R_obs, ...
                                           u_obs, v_obs, w_obs, ...
                                           LAMBDA_P, LAMBDA_S, ...
                                           LAMBDA_PRF, LAMBDA_OB, ...
                                           LAMBDA_T, LAMBDA_H, LAMBDA_F, ...
                                           t_array, t_diff, SAFE_DISTANCE, ...
                                           G, RHO, e, S, b, mass, AR, k, C_D_0, ...
                                           T_max, T_min, gamma_max, gamma_min, ...
                                           phi_max, phi_min, ...
                                           3);

            nonlcon = @(x) problem_constraints(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                               x_obs, y_obs, z_obs, R_obs, ...
                                               u_obs, v_obs, w_obs, ...
                                               R, R_dot, R_dot_dot, R_int, ...
                                               t_array, t_diff, NUM_LOCAL_WAYPOINTS, num_obs, ...
                                               x_d, y_d, z_d, psi_d, gamma_d, ...
                                               G, RHO, e, S, b, T_max, mass, AR, k, C_D_0, T_min, ...
                                               SAFE_DISTANCE, 3);
        else
            optimalWayPoints = [C_u(3:7); C_v(3:7); C_v(3:7)];
            objective = @(x) cost_function(x, C_u, C_v, C_w, ...
                                           x_0, y_0, z_0, ...
                                           R, R_dot, R_dot_dot, R_int, ...
                                           x_d, y_d, z_d, ...
                                           u_d, v_d, w_d, ...
                                           psi_d, gamma_d, ...
                                           HORIZON_TIME, NUM_LOCAL_WAYPOINTS, ...
                                           num_obs, x_obs, y_obs, z_obs, R_obs, ...
                                           u_obs, v_obs, w_obs, ...
                                           LAMBDA_P, LAMBDA_S, ...
                                           LAMBDA_PRF, LAMBDA_OB, ...
                                           LAMBDA_T, LAMBDA_H, LAMBDA_F, ...
                                           t_array, t_diff, SAFE_DISTANCE, ...
                                           G, RHO, e, S, b, mass, AR, k, C_D_0, ...
                                           T_max, T_min, gamma_max, gamma_min, ...
                                           phi_max, phi_min, ...
                                           2);

            nonlcon = @(x) problem_constraints(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                               x_obs, y_obs, z_obs, R_obs, ...
                                               u_obs, v_obs, w_obs, ...
                                               R, R_dot, R_dot_dot, R_int, ...
                                               t_array, t_diff, NUM_LOCAL_WAYPOINTS, num_obs, ...
                                               x_d, y_d, z_d, psi_d, gamma_d, ...
                                               G, RHO, e, S, b, T_max, mass, AR, k, C_D_0, T_min, ...
                                               SAFE_DISTANCE, 2);
        end

        %% Plot tajectory and prepare for next step
        [x, y, z, u, v, w, ...
         u_dot, v_dot, w_dot, ...
         u_dot_dot, v_dot_dot, w_dot_dot, ...
         V, V_dot, ...
         gamma, psi, psi_dot, gamma_dot, ...
         psi_dot_dot, gamma_dot_dot, ...
         phi, load_factor, C_L, C_D, D, T] = calc_traj_states(x_0, y_0, z_0, C_u, C_v, C_w, ...
                                                              R_int, R, R_dot, R_dot_dot, ...
                                                              t_array, t_diff, HORIZON_TIME, ...
                                                              G, RHO, S, C_D_0, k, mass);
        %{
        plot_trajectory(x, y, z, u, v, w, ...
                        u_dot, v_dot, w_dot, ...
                        u_dot_dot, v_dot_dot, w_dot_dot, ...
                        V, V_dot, gamma, psi, ...
                        psi_dot, psi_dot_dot, ...
                        gamma_dot, gamma_dot_dot, ...
                        phi, load_factor, T, t_array, ...
                        past_x, past_y, past_z, ...
                        num_obs, x_obs, y_obs, z_obs, R_obs);
        %}
        if any(x(1) >= x_obs(end) + R_obs(end)*20)
            break
        end

        % Define the initial conditions for the next step                                                 
        u_0 = u(NUM_TRAV_WAYPOINTS + 1);
        v_0 = v(NUM_TRAV_WAYPOINTS + 1);
        w_0 = w(NUM_TRAV_WAYPOINTS + 1);

        u_0_dot = u_dot(NUM_TRAV_WAYPOINTS + 1);
        v_0_dot = v_dot(NUM_TRAV_WAYPOINTS + 1);
        w_0_dot = w_dot(NUM_TRAV_WAYPOINTS + 1);

        u_0_dot_dot = u_dot_dot(NUM_TRAV_WAYPOINTS + 1);
        v_0_dot_dot = v_dot_dot(NUM_TRAV_WAYPOINTS + 1);
        w_0_dot_dot = w_dot_dot(NUM_TRAV_WAYPOINTS + 1);

        x_0 = x(NUM_TRAV_WAYPOINTS + 1);
        y_0 = y(NUM_TRAV_WAYPOINTS + 1);
        z_0 = z(NUM_TRAV_WAYPOINTS + 1);

        phi_current = phi(NUM_TRAV_WAYPOINTS + 1);
        T_current = T(NUM_TRAV_WAYPOINTS + 1);
        n_current = load_factor(NUM_TRAV_WAYPOINTS + 1);

        % Save initial points in past trajectory variables
        past_V = [past_V V(1:NUM_TRAV_WAYPOINTS)];
        past_x = [past_x x(1:NUM_TRAV_WAYPOINTS)];
        past_y = [past_y y(1:NUM_TRAV_WAYPOINTS)];
        past_z = [past_z z(1:NUM_TRAV_WAYPOINTS)];
        past_psi = [past_psi psi(1:NUM_TRAV_WAYPOINTS)];
        past_gamma = [past_gamma gamma(1:NUM_TRAV_WAYPOINTS)];
        past_u = [past_u u(1:NUM_TRAV_WAYPOINTS)];
        past_v = [past_v v(1:NUM_TRAV_WAYPOINTS)];
        past_w = [past_w w(1:NUM_TRAV_WAYPOINTS)];

        % Global trajectory
        % Remove the first points from the global trajectory
        V_g(1:NUM_TRAV_WAYPOINTS) = [];
        x_g(1:NUM_TRAV_WAYPOINTS) = [];
        y_g(1:NUM_TRAV_WAYPOINTS) = [];
        z_g(1:NUM_TRAV_WAYPOINTS) = [];
        psi_g(1:NUM_TRAV_WAYPOINTS) = [];
        gamma_g(1:NUM_TRAV_WAYPOINTS) = [];
        u_g(1:NUM_TRAV_WAYPOINTS) = [];
        v_g(1:NUM_TRAV_WAYPOINTS) = [];
        w_g(1:NUM_TRAV_WAYPOINTS) = [];

        if length(V_g) < NUM_LOCAL_WAYPOINTS
            u_g(1:NUM_LOCAL_WAYPOINTS) = u_g(end);
            v_g(1:NUM_LOCAL_WAYPOINTS) = v_g(end);
            w_g(1:NUM_LOCAL_WAYPOINTS) = w_g(end);
            V_g(1:NUM_LOCAL_WAYPOINTS) = V_g(end);
            x_g(1:NUM_LOCAL_WAYPOINTS) = linspace(x_g(1), x_g(1) + HORIZON_TIME*u_g(1), NUM_LOCAL_WAYPOINTS);
            y_g(1:NUM_LOCAL_WAYPOINTS) = linspace(y_g(1), y_g(1) + HORIZON_TIME*v_g(1), NUM_LOCAL_WAYPOINTS);
            z_g(1:NUM_LOCAL_WAYPOINTS) = linspace(z_g(1), z_g(1) + HORIZON_TIME*w_g(1), NUM_LOCAL_WAYPOINTS);
            psi_g(1:NUM_LOCAL_WAYPOINTS) = psi_g(end);
            gamma_g(1:NUM_LOCAL_WAYPOINTS) = gamma_g(end);
        end 

        % Update Demanded Trajectory
        V_d(1:NUM_LOCAL_WAYPOINTS) = V_g(1:NUM_LOCAL_WAYPOINTS);
        x_d(1:NUM_LOCAL_WAYPOINTS) = x_g(1:NUM_LOCAL_WAYPOINTS);
        y_d(1:NUM_LOCAL_WAYPOINTS) = y_g(1:NUM_LOCAL_WAYPOINTS);
        z_d(1:NUM_LOCAL_WAYPOINTS) = z_g(1:NUM_LOCAL_WAYPOINTS);
        psi_d(1:NUM_LOCAL_WAYPOINTS) = psi_g(1:NUM_LOCAL_WAYPOINTS);
        gamma_d(1:NUM_LOCAL_WAYPOINTS) = gamma_g(1:NUM_LOCAL_WAYPOINTS);
        u_d = u_g(1:NUM_LOCAL_WAYPOINTS);
        v_d = v_g(1:NUM_LOCAL_WAYPOINTS);
        w_d = w_g(1:NUM_LOCAL_WAYPOINTS);

        % Update obstacle position
        x_obs = x_obs + u_obs*(t_diff(1)*NUM_TRAV_WAYPOINTS);
        y_obs = y_obs + v_obs*(t_diff(1)*NUM_TRAV_WAYPOINTS);
        z_obs = z_obs + w_obs*(t_diff(1)*NUM_TRAV_WAYPOINTS);

        % toc
        % pause(0.5)
    end
    [energy, convergence, smoothness] = evaluate_hypercube_response(V, T, x, y, z, x_d, y_d, z_d, t_array);
    fprintf("For LAMBDA_P = %d, LAMBDA_S = %d, LAMBDA_PRF = %d and LAMBDA_T = %d, the energy was %f. Did the trajectory converged into the global trajectory? %s \n", LAMBDA_P, LAMBDA_S, LAMBDA_PRF, LAMBDA_T, energy, mat2str(convergence))
    if row == 1
        writematrix([LAMBDA_P LAMBDA_S LAMBDA_PRF LAMBDA_T energy convergence], 'response.xls')
    else
        writematrix([LAMBDA_P LAMBDA_S LAMBDA_PRF LAMBDA_T energy convergence], 'response.xls', 'WriteMode', 'append')
    end
end