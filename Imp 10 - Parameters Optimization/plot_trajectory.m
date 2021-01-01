function plot_trajectory(x, y, z, u, v, w, ...
                         u_dot, v_dot, w_dot, ...
                         u_dot_dot, v_dot_dot, w_dot_dot, ...
                         V, V_dot, gamma, psi, ...
                         psi_dot, psi_dot_dot, ...
                         gamma_dot, gamma_dot_dot, ...
                         phi, n, T, t_array, ...
                         past_x, past_y, past_z, ...
                         x_d, y_d, z_d, ...
                         m, x_obs, y_obs, z_obs, R_obs)
    
	% Plot position, acceleration, jerk, angles, load factor, drag and
    % thrust on each step of the trajectory                                            
    
    %% Two Dimensional Position
    %%{
    figure(1)
    subplot(3,1,1);
    plot(x, y, 'color','k','marker','.','markersize',16)
    for i = 1 : m
        viscircles([x_obs(i) y_obs(i)], R_obs(i));
    end
    grid on
    axis equal
    xlabel("x")
    ylabel("y")

    subplot(3,1,2);
    plot(x, z, 'color','k','marker','.','markersize',16)
    for i = 1 : m
        viscircles([x_obs(i) z_obs(i)], R_obs(i));
    end
    grid on
    axis equal
    xlabel("x")
    ylabel("z")
    
    subplot(3,1,3);
    plot(y, z, 'color','k','marker','.','markersize',16)
    for i = 1 : m
        viscircles([y_obs(i) z_obs(i)], R_obs(i));
    end
    grid on
    axis equal
    xlabel("y")
    ylabel("z")
    
    %%}
    
    %% Two Dimensional velocities
    
    %{
    figure(4)
    subplot(3,1,1);
    plot(t_array, u, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("u")

    %figure(5)
    subplot(3,1,2);
    plot(t_array, v, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("v")

    %figure(6)
    subplot(3,1,3);
    plot(t_array, w, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("w")

    figure(7)
    plot(t_array, V, 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("V")
    
    %}
    
    %% Two Dimensional acceleration

    %{
    figure(8)
    subplot(3,1,1);
    plot(t_array, u_dot, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("$\dot{u}$",'interpreter','latex')

    %figure(9)
    subplot(3,1,2);
    plot(t_array, v_dot, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("$\dot{v}$",'interpreter','latex')

    %figure(10)
    subplot(3,1,3);
    plot(t_array, w_dot, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("$\dot{w}$",'interpreter','latex')

   
    figure(11)
    plot(t_array, sqrt(u_dot.^2 + v_dot.^2 + w_dot.^2), 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("$\dot{V}$",'interpreter','latex')
    
    %figure(62)
    %plot(t_array(2:end), V_dot, 'color','k','marker','.','markersize',16)
    %xlabel("t")
    %ylabel("$\dot{V}$",'interpreter','latex')
    %}
 
 
    % Convert color code to 1-by-3 RGB array (0~1 each)
    str = '#D9FFFF';
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

    str = '#fbd9d3';
    prev_color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    
    figure(12)
    plot3(x, y, z,'-o','Color','b','MarkerSize',10, 'MarkerFaceColor', color)
    hold on
    plot3(past_x, past_y, past_z,'-o','Color','r','MarkerSize',10, 'MarkerFaceColor', prev_color)
    plot3(x_d, y_d, z_d,'-o','Color', [100, 100, 100]/255,'MarkerSize',10, 'MarkerFaceColor', [200, 200, 200]/255)
    plot_obstacle(m, x_obs, y_obs, z_obs, R_obs)
    grid on
    axis equal
    xlabel("x")
    ylabel("y")
    zlabel("z")
    %view(-101.9, 1.2) % Different point of view
    hold off

    %{
    figure(13)
    plot(t_array, z, 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("z")
    %}
    
    %{
    figure(14)
    %plot(t_array, [1 1 n], 'color','k','marker','.','markersize',16)
    plot(t_array, [1 n], 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("n")
    %}
    
    %%{
    figure(15)
    plot(t_array(2:end), T, 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("T")
    %%}