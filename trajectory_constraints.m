function [c, ceq] = trajectory_constraints(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                           x_obs, y_obs, z_obs, R_obs, ...
                                           u_obs, v_obs, w_obs, R_int, ...
                                           t_array, t_diff, n, m, ...
                                           safe_distance)
    %{
      Performs the calculation of nonlcon, which are the non-linear
      constraints, as defined in the fmincon function.
    %}
        
    C_u_ = [C_u(1:3); x(1:4)];
    C_v_ = [C_v(1:3); x(5:8)];
    C_w_ = [C_w(1:3); x(9:12)];  

    %% Calculate position
    x = x_0 + C_u_'*R_int;
    y = y_0 + C_v_'*R_int;
    z = z_0 + C_w_'*R_int;
    
    if m == 0 
        % No obstacles scenario
        c = 0;
    else
        % Iterate through all obstacles setting the constraints
        for i = 1 : m
            % Project future positions of the current obstacle, 
            % based on its current velocity
            x_cur_obs = x_obs(i) + u_obs(i)*t_array;
            y_cur_obs = y_obs(i) + v_obs(i)*t_array;
            z_cur_obs = z_obs(i) + w_obs(i)*t_array;

            % Calculate distance between the vehicle and the obstacle
            % at each point in time
            distance = sqrt((x_cur_obs - x).^2 + (y_cur_obs- y).^2 + (z_cur_obs- z).^2);

            c(1 + n*(i-1) : n*i) = R_obs(i) + safe_distance - distance;
        end
    end
   
    ceq = 0;
