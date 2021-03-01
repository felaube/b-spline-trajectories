function [convergence_cost, smoothness] = evaluate_hypercube_response(x, y, z, x_d, y_d, z_d, past_x, past_y, past_z)
    
    %{
      Calculate convergence cost and smoothness 
      relative to the trajectory
    %}
    convergence_cost = sum((x_d - x).^2 + (y_d - y).^2 + (z_d - z).^2);
    
    diff_x = diff(past_x);
    diff_y = diff(past_y);
    diff_z = diff(past_z);
    
    segments = [diff_x', diff_y', diff_z'];
    alpha = zeros(1, size(segments, 1) -1);
    
    for i = 1 : length(alpha)
        alpha(i) = acos(dot(-segments(i, :), segments(i+1, :))./(norm(-segments(i, :)).*norm(segments(i+1, :))));
    end
    
    smoothness = sum(alpha);
    