function [c, ceq] = trajectory_constraints(x, C_u, C_v, C_w, x_0, y_0, z_0, ...
                                           x_obs, y_obs, z_obs, R_obs, ...
                                           u_obs, v_obs, w_obs, ...
                                           R, R_dot, R_dot_dot, R_int, ...
                                           t_array, t_diff, n, m, ...
                                           x_d, y_d, z_d, psi_d, gamma_d, ...
                                           g, rho, e, S, b, T_max, mass, AR, k, C_D_0, T_min, ...
                                           safe_distance, num_set)
    
    % Performs the calculation of nonlcon, which are the non-linear
    % constraints, as defined for the fmincon function.

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

    %% Calculate position
    x = x_0 + C_u_'*R_int;
    y = y_0 + C_v_'*R_int;
    z = z_0 + C_w_'*R_int;
    
    %{
    %% Calculate speeds and their corresponding derivatives
    u = C_u_'*R;

    v = C_v_'*R;

    w = C_w_'*R;

    V = sqrt(u.^2 + v.^2 + w.^2);

    %% Calculate flight path angle, heading angle and their derivatives
    gamma = asin(w./V);
    psi = asin(v./(V.*cos(gamma)));

    gamma(isnan(gamma)) = 0; % Not sure if this is necessary. CHECK LATER
    psi(isnan(psi)) = 0; % Not sure if this is necessary. CHECK LATER

    
    % There are cases in which V_a is composed only by v_a.
    % In those cases, the value of v_a./(V_a.*cos(gamma_a)) can get
    % really close to 1, but not equal to 1. This can result in values
    % of psi with imaginary part. When that happens, we have to 
    % remove the imaginary part
    if any(imag(psi))
        psi = real(psi);
    end
    
    psi_dot = diff(psi)./t_diff;

    gamma_dot = diff(gamma)./t_diff;
    gamma_dot_dot = diff(gamma_dot)./t_diff(2:end);

    %% Calculate bank angle
    phi = atan2(psi_dot(1:end).*V(2:end).*cos(gamma(2:end)), g*cos(gamma(2:end)) + V(2:end).*gamma_dot);
    %}
    
    if m == 0
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
    
    %c = 0;
    ceq = 0;
    
    %ceq(1) = gamma(end) - gamma_d(end);
    %ceq(2) = psi(end) - psi_d(end);
    %ceq(3) = phi(end);
    %ceq(4) = x(end) - x_d(end);
    %ceq(5) = y(end) - y_d(end);
    %ceq(6) = z(end) - z_d(end);
