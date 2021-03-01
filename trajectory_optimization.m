function trajectory_evaluation = trajectory_optimization(lambdas, opts)
    
    if nargin < 2
        % Set default optimization settings
        opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
    end
    
    if nargin < 1
        % Set default lambdas values
        lambdas = [9, 10, 7, 6];
        
    end
        
    LAMBDA_P = lambdas(1); % position
    LAMBDA_S = lambdas(2); % speed
    LAMBDA_LIM = lambdas(3); % vehicle contraints (limits)
    LAMBDA_AT = lambdas(4); % terminal attitude

    % Discretization of the receding horizon trajectory

    % Number of waypoints in local trajectory
    NUM_LOCAL_WAYPOINTS = 51;
    % Number of traveled points before running the optimization process again
    NUM_TRAV_WAYPOINTS = 4;
    % Receding horizon time
    HORIZON_TIME = 100;

    tal = linspace(0, 1, NUM_LOCAL_WAYPOINTS);
    t_array = tal*HORIZON_TIME;

    % t_diff, as implemented here, is an array. It could just as well
    % be a constant, since the points are "equitemporal". However,
    % by maintaining it as an array, the script is ready to deal
    % with non "equitemporal" points, if the necessity arises
    t_diff = diff(t_array);

    % Obstacles

    % The method does not behave well when the starting position 
    % of the obstacle is too close to the starting position of the vehicle
    OBSTACLES_SCENARIO_ID = 8.1;
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
    % In all scenarios of this implementation, the global trajectory
    % is the same, and it is a straight level flight, so it is defined 
    % with only NUM_LOCAL_WAYPOINTS points, but it could be more
    V_g(1:NUM_LOCAL_WAYPOINTS) = 30;
    psi_g(1:NUM_LOCAL_WAYPOINTS) = 0;
    gamma_g(1:NUM_LOCAL_WAYPOINTS) = 0;
    
    u_g = V_g.*cos(gamma_g).*cos(psi_g);
    v_g = V_g.*cos(gamma_g).*sin(psi_g);
    w_g = V_g.*sin(gamma_g);

    x_g = 0 : t_diff(1) * u_g(1) : HORIZON_TIME * u_g(1);
    y_g = 1000 * ones(1, length(x_g));
    z_g = 1000 * ones(1, length(x_g));
    
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
    y_0 = 1000;
    z_0 = 1200;

    % Vector of Coefficients
    C_u = zeros(7, 1);
    C_v = zeros(7, 1);
    C_w = zeros(7, 1);

    C_u(4:7) = 30;
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
    past_x_d = [];
    past_y_d = [];
    past_z_d = [];
    past_psi = [];
    past_phi = [];
    past_gamma = [];
    past_u = [];
    past_v = [];
    past_w = [];
    past_T = [];
    
    % Initialize variables used to store the objective function value
    % and the optmization results
    optim_vars = cell(3, 28);
    cost_val(1:28) = Inf;
    
    while true

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

        initial_condition = [C_u(4:7); C_v(4:7); C_w(4:7)];

        objective = @(x) cost_function(x, C_u, C_v, C_w, ...
                                       x_0, y_0, z_0, ...
                                       R, R_dot, R_dot_dot, R_int, ...
                                       x_d, y_d, z_d, ...
                                       u_d, v_d, w_d, ...
                                       psi_d, gamma_d, ...
                                       HORIZON_TIME, NUM_LOCAL_WAYPOINTS, ...
                                       R_obs, ...
                                       LAMBDA_P, LAMBDA_S, ...
                                       LAMBDA_LIM, LAMBDA_AT, ...
                                       t_array, t_diff, SAFE_DISTANCE, ...
                                       G, RHO, e, S, b, mass, AR, k, C_D_0, ...
                                       T_max, T_min, gamma_max, gamma_min, ...
                                       phi_max, phi_min);

        nonlcon = @(x) trajectory_constraints(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                              x_obs, y_obs, z_obs, R_obs, ...
                                              u_obs, v_obs, w_obs, R_int, ...
                                              t_array, t_diff, NUM_LOCAL_WAYPOINTS, num_obs, ...
                                              SAFE_DISTANCE);

        [optimalWayPoints, fval, ~, output] = fmincon(objective, initial_condition(:), [],[],[],[],[],[],nonlcon,opts);

        C_u = [C_u(1); C_u(2); C_u(3); optimalWayPoints(1:4)];
        C_v = [C_v(1); C_v(2); C_v(3); optimalWayPoints(5:8)];
        C_w = [C_w(1); C_w(2); C_w(3); optimalWayPoints(9:12)];
        
        %% Plot tajectory and prepare for next step
        [x, y, z, u, v, w, ...
         u_dot, v_dot, w_dot, ...
         u_dot_dot, v_dot_dot, w_dot_dot, ...
         V, V_dot, ...
         gamma, psi, psi_dot, gamma_dot, ...
         psi_dot_dot, gamma_dot_dot, ...
         phi, load_factor, C_L, C_D, D, T] = calculate_trajectory_states(x_0, y_0, z_0, C_u, C_v, C_w, ...
                                                                         R_int, R, R_dot, R_dot_dot, ...
                                                                         t_array, t_diff, HORIZON_TIME, ...
                                                                         G, RHO, S, C_D_0, k, mass);
        
        plot_trajectory(x, y, z, u, v, w, ...
                        u_dot, v_dot, w_dot, ...
                        u_dot_dot, v_dot_dot, w_dot_dot, ...
                        V, V_dot, gamma, psi, ...
                        psi_dot, psi_dot_dot, ...
                        gamma_dot, gamma_dot_dot, ...
                        phi, load_factor, T, t_array, ...
                        past_x, past_y, past_z, ...
                        past_x_d, past_y_d, past_z_d, ...
                        past_u, past_v, past_w, ...
                        past_phi, past_gamma, past_T, ...
                        x_d, y_d, z_d, ...
                        num_obs, x_obs, y_obs, z_obs, R_obs);
        
        % Uncomment line below to add a delay to the optimization
        % pause(0.25)
        
        % Check trajectory stop condition
        if R_obs(end) == 0 % Stop condition for obstacle free scenario
            if x(end) >= 7000
                break
            end
        else
            if x(1) >= (x_obs(end) + R_obs(end)*20) || length(x) > 2000
                break
            end
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

        % Save initial points in past trajectory variables
        past_V = [past_V V(1:NUM_TRAV_WAYPOINTS)];
        past_x = [past_x x(1:NUM_TRAV_WAYPOINTS)];
        past_y = [past_y y(1:NUM_TRAV_WAYPOINTS)];
        past_z = [past_z z(1:NUM_TRAV_WAYPOINTS)];
        past_x_d = [past_x_d x_d(1:NUM_TRAV_WAYPOINTS)];
        past_y_d = [past_y_d y_d(1:NUM_TRAV_WAYPOINTS)];
        past_z_d = [past_z_d z_d(1:NUM_TRAV_WAYPOINTS)];
        past_phi = [past_phi phi(1:NUM_TRAV_WAYPOINTS)];
        past_psi = [past_psi psi(1:NUM_TRAV_WAYPOINTS)];
        past_gamma = [past_gamma gamma(1:NUM_TRAV_WAYPOINTS)];
        past_u = [past_u u(1:NUM_TRAV_WAYPOINTS)];
        past_v = [past_v v(1:NUM_TRAV_WAYPOINTS)];
        past_w = [past_w w(1:NUM_TRAV_WAYPOINTS)];
        past_T = [past_T T(1:NUM_TRAV_WAYPOINTS)];

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

        if length(x_g) < NUM_LOCAL_WAYPOINTS
            % Extend global trajectory, considering that the current
            % states should not change
            x_g(1:NUM_LOCAL_WAYPOINTS) = linspace(x_g(1), x_g(1) + HORIZON_TIME*u_g(1), NUM_LOCAL_WAYPOINTS);
            y_g(1:NUM_LOCAL_WAYPOINTS) = linspace(y_g(1), y_g(1) + HORIZON_TIME*v_g(1), NUM_LOCAL_WAYPOINTS);
            z_g(1:NUM_LOCAL_WAYPOINTS) = linspace(z_g(1), z_g(1) + HORIZON_TIME*w_g(1), NUM_LOCAL_WAYPOINTS);
            u_g(1:NUM_LOCAL_WAYPOINTS) = u_g(end);
            v_g(1:NUM_LOCAL_WAYPOINTS) = v_g(end);
            w_g(1:NUM_LOCAL_WAYPOINTS) = w_g(end);
            V_g(1:NUM_LOCAL_WAYPOINTS) = V_g(end);
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
        x_obs = x_obs + u_obs*(t_array(NUM_TRAV_WAYPOINTS + 1));
        y_obs = y_obs + v_obs*(t_array(NUM_TRAV_WAYPOINTS + 1));
        z_obs = z_obs + w_obs*(t_array(NUM_TRAV_WAYPOINTS + 1));
        
    end

    [convergence, smoothness] = evaluate_hypercube_response(x, y, z, x_d, y_d, z_d, past_x, past_y, past_z);
    fprintf("For LAMBDA_P = %d, LAMBDA_S = %d, LAMBDA_PRF = %d and LAMBDA_T = %d. Smoothness: %f. Convergence cost: %f \n", LAMBDA_P, LAMBDA_S, LAMBDA_LIM, LAMBDA_AT, smoothness, convergence)
  
    trajectory_evaluation = [convergence, smoothness];

