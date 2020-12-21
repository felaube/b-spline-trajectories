function [energy, convergence, smoothness] = evaluate_hypercube_response(V, T, x, y, z, x_d, y_d, z_d, t_array)

    energy = trapz(t_array(2:end), T.*V(2:end));
    distance_to_required = sum(abs(x - x_d)./x_d + abs(y - y_d)./y_d + abs(z - z_d)./z_d);
    
    convergence = (distance_to_required <= (0.05*length(x)));
    
    smoothness = 0;