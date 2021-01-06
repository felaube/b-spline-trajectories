function [energy, convergence, smoothness] = evaluate_hypercube_response(V, T, x, y, z, x_d, y_d, z_d, past_x, past_y, past_z, t_array)

    energy = trapz(t_array(2:end), T.*V(2:end));
    % distance_to_required = zeros(1, length(x)*3);
    % distance_to_required(1:length(x)) = abs(x - x_d)./x_d;
    % distance_to_required(1+length(x):2*length(x)) = abs(y - y_d)./y_d;
    % distance_to_required(1+2*length(x):3*length(x)) = abs(z - z_d)./z_d;
    
    % convergence = all(distance_to_required <= 0.05);
    
    % convergence = sum(distance_to_required);
    
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
    % angle_x_y = asin(diff_y./diff_x);
    % angle_x_z = asin(diff_z./diff_x);
    % angle_y_z = asin(diff_z./diff_y);
    
    % diff_angle_x_y = diff(angle_x_y);
    % diff_angle_x_z = diff(angle_x_z);
    % diff_angle_y_z = diff(angle_y_z);
    
    % smoothness = sum(abs(diff_angle_x_y) + abs(diff_angle_x_z)+ abs(diff_angle_y_z));