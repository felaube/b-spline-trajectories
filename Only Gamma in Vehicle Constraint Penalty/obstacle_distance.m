function distance = obstacle_distance(x_obs, y_obs, z_obs, R_obs, x_a, y_a, z_a)
    
    % Find the trajectory closest point to the obstacle
    
    % Get a generic sphere
    % [X,Y,Z] = sphere;
    
    % Update the points to match the demanded obstacle
    % X = x_obs + X * R_obs;
    % Y = y_obs + Y * R_obs;
    % Z = z_obs + Z * R_obs;
    
    safe_distance = 50;
    
    possible_distances = sqrt((x_obs - x_a).^2 + (y_obs - y_a).^2 + (z_obs - z_a).^2) - R_obs - safe_distance;
    
    possible_distances(possible_distances<=0) = 1e-8;
    
    distance = min(possible_distances);
    
    %distance = 1e100;
    