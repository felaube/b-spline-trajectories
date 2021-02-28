function performance_margin = vehicle_constraints(V, V_dot, ...
                                                  gamma, ...
                                                  gamma_dot, psi_dot, ...
                                                  gamma_dot_dot)
    % Aircraft: DECODE 2
    
                                              
    %% Constants 
    % Oswald factor [-]
    % e = 0.8;
    
    % Wing area [m²]
    S = 1.554;
    
    % Span [m]
    b = 3.9422;
    
    % Maximum thrust (all operating engines) [N]
    T_max = 130;
    
    % Aircraft mass [kg]
    m = 25;
    
    % Gravitational acceleration [m/s²]
    g = 9.81;
    
    % Air density (speed is considered to be AES) [kg/m³]
    rho = 1.225;
    
    % Aspect Ratio [-]
    %A = b^2/S;
    
    % k value [-]
    k = 0.0334;
    
    % Cd0 value (guessed) [-]
    C_D_0 = 0.0715;
    
    %% Calculate bank angle
    phi = atan2(psi_dot(2:end).*V(3:end).*cos(gamma(3:end)), g*cos(gamma(3:end)) + V(3:end).*gamma_dot_dot);
    
    %% Calculate load factor
    n = (g*cos(gamma(3:end)) + V(3:end).*gamma_dot(2:end))./(g*cos(phi));
    
    %% Calculate thrust
    C_L = 2*n*m*g./(rho*S*V(3:end).^2);
    C_D = C_D_0 + k*C_L.^2;
    
    D = 1/2*rho*S*C_D.*V(3:end).^2;
    
    T = D + m*V_dot(2:end) + m*g*sin(gamma(3:end));
    
    performance_margin = 100 - 100*(max(T)/T_max);
    
    performance_margin = 1e100;
    
    % performance_margin = 100 - 100*(max(gamma)/deg2rad(20));