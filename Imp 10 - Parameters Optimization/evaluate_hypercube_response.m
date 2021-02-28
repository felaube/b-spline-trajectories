function [energy, convergence, smoothness, negative_T] = evaluate_hypercube_response(V, T, x, y, z, x_d, y_d, z_d, past_x, past_y, past_z, t_array, past_T)

    energy = trapz(t_array(2:end), T.*V(2:end));
    
    convergence = sum((x_d - x).^2 + (y_d - y).^2 + (z_d - z).^2);
    
    diff_x = diff(past_x);
    diff_y = diff(past_y);
    diff_z = diff(past_z);
    
    segments = [diff_x', diff_y', diff_z'];
    alpha = zeros(1, size(segments, 1) -1);
    
    for i = 1 : length(alpha)
        alpha(i) = acos(dot(-segments(i, :), segments(i+1, :))./(norm(-segments(i, :)).*norm(segments(i+1, :))));
    end
    
    smoothness = sum(alpha);
    
    negative_T = any([T past_T] < 0);
    