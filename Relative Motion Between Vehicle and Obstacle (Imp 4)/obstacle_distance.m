function distance = obstacle_distance(x_obs, y_obs, z_obs, R_obs, ...
                                      x_a, y_a, z_a, safe_distance)

    %{
        Calculate the distance from the obstacle to each waypoint of the
        trajectory and returns the minimum distance. I.e. returns the 
        distance between the obstacle and the closest trajectory point to it.
    %}


    possible_distances = sqrt((x_obs - x_a).^2 + (y_obs - y_a).^2 + (z_obs - z_a).^2) - R_obs - safe_distance;

    possible_distances(possible_distances<=0) = 1e-8;

    distance = min(possible_distances);
    