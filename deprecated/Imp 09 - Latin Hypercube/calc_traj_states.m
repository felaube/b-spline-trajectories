function [x, y, z, u, v, w, ...
          u_dot, v_dot, w_dot, ...
          u_dot_dot, v_dot_dot, w_dot_dot, ...
          V, V_dot, ...
          gamma, psi, psi_dot, gamma_dot, ...
          psi_dot_dot, gamma_dot_dot, ...
          phi, n, C_L, C_D, D, T] = calc_traj_states(x_0, y_0, z_0, C_u, C_v, C_w, ...
                                                     R_int, R, R_dot, R_dot_dot, ...
                                                     t_array, t_diff, t_h, ...
                                                     g, rho, S, C_D_0, k, mass)
                                                 
    % Calculate position, acceleration, jerk, angles, load factor, drag and
    % thrust on each step of the trajectory
    
    %% Calculate position
    x = x_0 + C_u'*R_int;
    y = y_0 + C_v'*R_int;
    z = z_0 + C_w'*R_int;
    
    %% Calculate velocities and their corresponding derivatives
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
    % V_dot = diff(V)./t_diff;
    V_dot = sqrt(u_dot.^2 + v_dot.^2 + w_dot.^2);
    %% Calculate flight path angle, heading angle and their derivatives
    gamma = asin(w./V);
    gamma(isnan(gamma)) = 0; % Not sure if this is necessary. CHECK LATER
    psi = asin(v./(V.*cos(gamma)));
    psi(isnan(psi)) = 0; % Not sure if this is necessary. CHECK LATER
    
    % There are cases in which V is composed only by v.
    % In those cases, the value of v./(V.*cos(gamma)) can get
    % really close to 1, but not equal to 1. This can result in values
    % of psi with imaginary part. When that happens, we have to 
    % remove the imaginary part
    if any(imag(psi))
        psi = real(psi);
    end
    
    psi_dot = diff(psi)./t_diff;
    psi_dot_dot = diff(psi_dot)./t_diff(2:end);

    gamma_dot = diff(gamma)./t_diff;
    gamma_dot_dot = diff(gamma_dot)./t_diff(2:end);


    %% Calculate bank angle
    phi = atan2(psi_dot(1:end).*V(2:end).*cos(gamma(2:end)), g*cos(gamma(2:end)) + V(2:end).*gamma_dot);

    %% Calculate load factor
    n = (g*cos(gamma(2:end)) + V(2:end).*gamma_dot(1:end))./(g*cos(phi));

    %% Calculate thrust
    C_L = 2*n*mass*g./(rho*S*V(2:end).^2);
    C_D = C_D_0 + k*C_L.^2;

    D = 1/2*rho*S*C_D.*V(2:end).^2;

    % T = D + mass*V_dot(1:end) + mass*g*sin(gamma(2:end));
    T = D + mass*V_dot(2:end) + mass*g*sin(gamma(2:end));