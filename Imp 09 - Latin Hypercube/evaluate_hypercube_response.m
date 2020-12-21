function [energy, convergence, smoothness] = evaluate_hypercube_response(V, T, x, y, z, x_d, y_d, z_d, t_array)

    energy = trapz(t_array(2:end), T.*V(2:end));
    distance_to_required = zeros(1, length(x)*3);
    distance_to_required(1:length(x)) = abs(x - x_d)./x_d;
    distance_to_required(1+length(x):2*length(x)) = abs(y - y_d)./y_d;
    distance_to_required(1+2*length(x):3*length(x)) = abs(z - z_d)./z_d;
    
    convergence = all(distance_to_required <= 0.05);
    
    smoothness = 0;