function trajectory_cost = cost_function(x, C_u, C_v, C_w, ...
                                         x_0, y_0, z_0, ...
                                         R, R_dot, R_dot_dot, R_int, ...
                                         x_d, y_d, z_d, ...
                                         u_d, v_d, w_d, ...
                                         psi_d, gamma_d, ...
                                         t_h, n, ...
                                         m, x_obs, y_obs, z_obs, R_obs, ...
                                         u_obs, v_obs, w_obs, ...
                                         lambda_p, lambda_s, ...
                                         lambda_prf, lambda_ob, ...
                                         lambda_t, lambda_h, lambda_f, ...
                                         t_array, t_diff, safe_distance, ...
                                         g, rho, e, S, b, mass, AR, k, C_D_0, ...
                                         T_max, T_min, gamma_max, gamma_min, ...
                                         phi_max, phi_min, ...
                                         num_set)

    % Calculate trajectory cost (objective function).
    
    if num_set == 3 
        % Configuration used in the first optimization process 
        C_u_ = [C_u(1:3); x(1:4)];
        C_v_ = [C_v(1:3); x(5:8)];
        C_w_ = [C_w(1:3); x(9:12)];
    else
        if num_set == 2
            % Configuration used in the optimization processes
            % that exist to avoid local minima.
            C_u_ = [C_u(1:2); x(1:5)];
            C_v_ = [C_v(1:2); x(6:10)];
            C_w_ = [C_w(1:2); x(11:15)];
        end
    end

	[x_a, y_a, z_a, u_a, v_a, w_a, ...
     u_a_dot, v_a_dot, w_a_dot, ...
     u_a_dot_dot, v_a_dot_dot, w_a_dot_dot, ...
     V_a, V_a_dot, ...
     gamma_a, psi_a, psi_a_dot, gamma_a_dot, ...
     psi_a_dot_dot, gamma_a_dot_dot, ...
     phi_a, n_a, C_L, C_D, D_a, T_a] = calculate_trajectory_states(x_0, y_0, z_0, C_u_, C_v_, C_w_, ...
                                                                   R_int, R, R_dot, R_dot_dot, ...
                                                                   t_array, t_diff, t_h, ...
                                                                   g, rho, S, C_D_0, k, mass);

    %% Calculate costs
    position_cost = sum((x_d - x_a).^2 + (y_d - y_a).^2 + (z_d - z_a).^2)/(n*((max(R_obs) + safe_distance)^2));
    speed_cost = sum((u_d - u_a).^2 + (v_d - v_a).^2 + (w_d - w_a).^2)/(n*(mean(u_d)^2 + mean(v_d)^2 + mean(w_d)^2));
    
    obstacles_constraints_penalty = 0;
    
    A_ob = 1;
    alpha_ob = 0.45;
    % Calculate distance to each obstacle and add to obstacles_constraints_penalty
    for i = 1 : m
        distance = obstacle_distance(x_obs(i), y_obs(i), z_obs(i), R_obs(i), ...
                                     u_obs(i), v_obs(i), w_obs(i), ...
                                     x_a, y_a, z_a, ...
                                     t_array, safe_distance);
        
        
        obstacles_constraints_penalty = obstacles_constraints_penalty + ...
                                                A_ob*exp(alpha_ob*(R_obs(i) + safe_distance - distance));
        
    end
    
    % Calculate vehicle constraints penalty
    %vehicle_constraints_penalty = (max(0, T_a - T_max) + ...
    %                               max(0, T_min - T_a) + ...
    %                               max(0, gamma_a(1:end-1) - gamma_max) + ...
    %                               max(0, gamma_min - gamma_a(1:end-1)) + ...
    %                               max(0, phi_a - phi_max) + ...
    %                               max(0, phi_min - phi_a)).^2;
                               
     vehicle_constraints_penalty = (max(0, T_a - T_max) + ...
                                    max(0, gamma_a(1:end-1) - gamma_max) + ...
                                    max(0, gamma_min - gamma_a(1:end-1)) + ...
                                    max(0, phi_a - phi_max) + ...
                                    max(0, phi_min - phi_a)).^2;
    
    vehicle_constraints_penalty = sum(vehicle_constraints_penalty);
    
    terminal_tracking_cost = lambda_h*((psi_d(n) - psi_a(n))/pi)^2 + lambda_f*((gamma_d(n) - gamma_a(n))/pi)^2;
    
    terminal_position_cost = 1e4*(1*(x_a(n) - x_d(n))^2 + 1*(y_a(n) - y_d(n))^2 + 1*(z_a(n) - z_d(n))^2);
    
    % Saturation of costs to avoid NaN and Inf results 
    % in the objective function
    position_cost = min(position_cost, 1e200);
    speed_cost = min(speed_cost, 1e200);
    vehicle_constraints_penalty = min(vehicle_constraints_penalty, 1e200);
    obstacles_constraints_penalty = min(obstacles_constraints_penalty, 1e200);
    terminal_tracking_cost = min(terminal_tracking_cost, 1e200);
    
    
    trajectory_cost = lambda_p*position_cost ... 
                    + lambda_s*speed_cost ...
                    + lambda_prf*vehicle_constraints_penalty ...
                    + lambda_ob*obstacles_constraints_penalty ...
                    + lambda_t*terminal_tracking_cost + 0*terminal_position_cost;