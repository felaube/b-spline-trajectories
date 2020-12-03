function performance_margin = vehicle_constraints_cost(gamma, gamma_dot, gamma_dot_dot, psi_dot, V, V_dot, ...
                                                       g, rho, e, S, b, mass, AR, k, C_D_0, ...
                                                       T_max, T_min, gamma_max, gamma_min);
    %{
        Calculate the performance margin, required for the 
        vehicle constraints penalty function, which is calculated
        in the objective function.
    %}
                                                       
    % Calculate bank angle
    phi = atan2(psi_dot(1:end).*V(2:end).*cos(gamma(2:end)), g*cos(gamma(2:end)) + V(2:end).*gamma_dot);

    % Calculate load factor
    n = (g*cos(gamma(2:end)) + V(2:end).*gamma_dot(1:end))./(g*cos(phi));

    % Calculate thrust
    C_L = 2*n*mass*g./(rho*S*V(2:end).^2);
    C_D = C_D_0 + k*C_L.^2;

    D = 1/2*rho*S*C_D.*V(2:end).^2;

    T = D + mass*V_dot(1:end) + mass*g*sin(gamma(2:end));

    T(isnan(T)) = 0; % Not sure if this is necessary. CHECK LATER

    % performance_margin = 1.5*(max(gamma)/gamma_max)^2 + 1.5*(max(T)/T_max);
    
    performance_margin = 100 - (T/T_max)*100;

    performance_margin_min = 1.5*(T_min/T_max);

    performance_margin(performance_margin<performance_margin_min) = performance_margin_min;