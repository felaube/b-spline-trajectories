function trajectory_cost = cost_function(x, C_u, C_v, C_w, ...
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
                                             t_array, t_diff)
    
    %% Construct vectors of coefficients 
    C_u_ = [C_u(1:3); x(1:4)];
    C_v_ = [C_v(1:3); x(5:8)];
    C_w_ = [C_w(1:3); x(9:12)];

    %% Calculate position
    x_a = x_0 + C_u_'*R_int;
    y_a = y_0 + C_v_'*R_int;
    z_a = z_0 + C_w_'*R_int;

    %% Calculate speeds and their corresponding derivatives
    u_a = C_u_'*R;
    u_a_dot = C_u_'*R_dot*1/t_h;
    u_a_dot_dot = C_u_'*R_dot_dot*(1/(t_h^2));

    v_a = C_v_'*R;
    v_a_dot = C_v_'*R_dot*1/t_h;
    v_a_dot_dot = C_v_'*R_dot_dot*(1/(t_h^2));

    w_a = C_w_'*R;
    w_a_dot = C_w_'*R_dot*1/t_h;
    w_a_dot_dot = C_w_'*R_dot_dot*(1/(t_h^2));

    V_a = sqrt(u_a.^2 + v_a.^2 + w_a.^2);
    V_a_dot = diff(V_a)./t_diff;
    
    %% Calculate flight path angle, heading angle and their derivatives
    gamma_a = asin(w_a./V_a);
    gamma_a(isnan(gamma_a)) = 0; % Not sure if this is necessary. CHECK LATER
    
    psi_a = asin(v_a./(V_a.*cos(gamma_a)));
    psi_a(isnan(psi_a)) = 0; % Not sure if this is necessary. CHECK LATER
    
    psi_a_dot = diff(psi_a)./t_diff;
    psi_a_dot_dot = diff(psi_a_dot)./t_diff(2:end);
    
    gamma_a_dot = diff(gamma_a)./t_diff;
    gamma_a_dot_dot = diff(gamma_a_dot)./t_diff(2:end);

    %% Calculate costs
    position_cost = sum((x_d - x_a).^2 + (y_d - y_a).^2 + (z_d - z_a).^2);
    speed_cost = sum((u_d - u_a).^2 + (v_d - v_a).^2 + (w_d - w_a).^2);
    vehicle_constraints_penalty = 0;
    obstacles_constraints_penalty = 0;
    
    A_ob = 1;
    % OBS: alpha_ob tem pouca influ�ncia no resultado entre 3 e 0.125,
    % 0.12 � OK, qualquer valor menor do que isso causa superposi��o
    % da aeronave e do obst�culo
    alpha_ob = 1;
    % Calculate distance to each obstacle and add to obstacles_constraints_penalty
    for i = 1 : m
        distance = obstacle_distance(x_obs(i), y_obs(i), z_obs(i), R_obs, x_a, y_a, z_a);
        
        %obstacles_constraints_penalty = obstacles_constraints_penalty + ...
        %                                    A_ob*1/distance;
        
        obstacles_constraints_penalty = obstacles_constraints_penalty + ...
                                                A_ob*exp(-alpha_ob*distance)/distance;
        
    end
    
    % Calculate vehicle constraints penalty
    
    A_p = 1;
    alpha_p = 1;
    
    performance_margin = vehicle_constraints(V_a, V_a_dot, ...
                                             gamma_a, ...
                                             gamma_a_dot, psi_a_dot, ...
                                             gamma_a_dot_dot);
    if distance < 0
        vehicle_constraints_penalty = A_p * exp(-alpha_p*performance_margin)/(-performance_margin);
    else
        vehicle_constraints_penalty = A_p * exp(-alpha_p*performance_margin)/performance_margin;
    end
    
    
    terminal_tracking_cost = lambda_h*(psi_d(n) - psi_a(n))^2 + lambda_f*(gamma_d(n) - gamma_a(n))^2;
    
        
    trajectory_cost = lambda_p*position_cost ... 
                    + lambda_s*speed_cost ...
                    + lambda_prf*vehicle_constraints_penalty ...
                    + lambda_ob*obstacles_constraints_penalty ...
                    + lambda_t*terminal_tracking_cost + 0*terminal_position_cost;