    function performance_margin = vehicle_constraints_cost(gamma, gamma_dot, gamma_dot_dot, gamma_min, gamma_max, ...
                                                       psi_dot, V, V_dot, ...
                                                       g, rho, e, S, b, mass, AR, k, C_D_0, ...
                                                       T, T_min, T_max, phi, n)
        
	% Calculate the performance margin, required for the 
    % vehicle constraints penalty function, which is calculated
    % in the objective function.

    % performance_margin = 1.5*(max(gamma)/gamma_max)^2 + 1.5*(max(T)/T_max);
    
    %{
    performance_margin = 100 - (T/T_max)*100;
    
    performance_margin_min = 100 - (T_min/T_max)*100;
    
    performance_margin_max = 90;

    performance_margin(performance_margin<performance_margin_min) = performance_margin_min;
    performance_margin(performance_margin>performance_margin_max) = performance_margin_max;
    %}
    
    performance_margin = 100 - (T/T_max)*100;
    