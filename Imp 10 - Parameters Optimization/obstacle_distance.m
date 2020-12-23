function distance = obstacle_distance(x_obs, y_obs, z_obs, R_obs, ...
                                      u_obs, v_obs, w_obs, ... 
                                      x_a, y_a, z_a, ...
                                      t_array, safe_distance)

    %{
        Calculate the distance from the obstacle to each waypoint of the
        trajectory and returns the minimum distance. I.e. returns the 
        distance between the obstacle and the closest trajectory point to it.
    %}

    % Project future positions of the obstacle, 
    % based on its current velocity
    x_obs = x_obs + u_obs*t_array;
    y_obs = y_obs + v_obs*t_array;
    z_obs = z_obs + w_obs*t_array;
    
    possible_distances = sqrt((x_obs - x_a).^2 + (y_obs - y_a).^2 + (z_obs - z_a).^2);

    distance = min(possible_distances);
    