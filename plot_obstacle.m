function plot_obstacle(m, x_obs, y_obs, z_obs, R_obs)
    
    %{
        Plot a sphere in the current figure
    %}
    [X,Y,Z] = sphere;
    
    for i = 1 : m
        hold on
        X2 = X * R_obs(i);
        Y2 = Y * R_obs(i);
        Z2 = Z * R_obs(i);
        surf(X2 + x_obs(i), Y2 + y_obs(i), Z2 + z_obs(i));        
        hold off
    end
