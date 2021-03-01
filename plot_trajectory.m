function plot_trajectory(x, y, z, u, v, w, ...
                         u_dot, v_dot, w_dot, ...
                         u_dot_dot, v_dot_dot, w_dot_dot, ...
                         V, V_dot, gamma, psi, ...
                         psi_dot, psi_dot_dot, ...
                         gamma_dot, gamma_dot_dot, ...
                         phi, n, T, t_array, ...
                         past_x, past_y, past_z, ...
                         past_x_d, past_y_d, past_z_d, ...
                         past_u, past_v, past_w, ...
                         past_phi, past_gamma, past_T,...
                         x_d, y_d, z_d, ...
                         m, x_obs, y_obs, z_obs, R_obs)
                     
	% Plot position, acceleration, jerk, angles, load factor, drag and
    % thrust on each step of the trajectory                                            
    
    fig_num = 0;
    %% Two Dimensional Position
    
    fig_num = fig_num + 1;
    figure(fig_num)
    subplot(2,1,1);
    plot(x, y, 'color','k','marker','.','markersize',16)
    for i = 1 : m
        viscircles([x_obs(i) y_obs(i)], R_obs(i));
    end
    grid on
    axis equal
    xlabel("x")
    ylabel("y")

    subplot(2,1,2);
    plot(x, z, 'color','k','marker','.','markersize',16)
    for i = 1 : m
        viscircles([x_obs(i) z_obs(i)], R_obs(i));
    end
    grid on
    axis equal
    xlabel("x")
    ylabel("z")
    
    %% Two Dimensional velocities
    
    
    fig_num = fig_num + 1;
    figure(fig_num)
    subplot(3,1,1);
    plot(t_array, u, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("u")

    subplot(3,1,2);
    plot(t_array, v, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("v")

    subplot(3,1,3);
    plot(t_array, w, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("w")

    fig_num = fig_num + 1;
    figure(fig_num)
    plot(t_array, V, 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("V")
    
    
    
    %% Two Dimensional acceleration

    
    fig_num = fig_num + 1;
    figure(fig_num)
    subplot(3,1,1);
    plot(t_array, u_dot, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("$\dot{u}$",'interpreter','latex')

    subplot(3,1,2);
    plot(t_array, v_dot, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("$\dot{v}$",'interpreter','latex')

    subplot(3,1,3);
    plot(t_array, w_dot, 'color','k','marker','.','markersize',16)
    grid on
    xlabel("t")
    ylabel("$\dot{w}$",'interpreter','latex')

   
    fig_num = fig_num + 1;
    figure(fig_num)
    plot(t_array, sqrt(u_dot.^2 + v_dot.^2 + w_dot.^2), 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("$\dot{V}$",'interpreter','latex')
    
    %%  Two Dimensional Historical Position
    
    fig_num = fig_num + 1;
    figure(fig_num)
    subplot(2,1,1);
    plot([past_x x], [past_y y], 'color','k','marker','.','markersize',16)
    for i = 1 : m
        viscircles([x_obs(i) y_obs(i)], R_obs(i));
    end
    set(gcf,'color','w');
    grid on
    axis equal
    xlabel("x [m]")
    ylabel("y [m]")

    subplot(2,1,2);
    plot([past_x x], [past_z z], 'color','k','marker','.','markersize',16)
    for i = 1 : m
        viscircles([x_obs(i) z_obs(i)], R_obs(i));
    end
    set(gcf,'color','w');
    grid on
    axis equal
    xlabel("x [m]")
    ylabel("z [m]")
    
    %%  Two Dimensional Historical Velocities
    
    fig_num = fig_num + 1;
    figure(fig_num)
    subplot(3,1,1);
    plot([past_x x], [past_u u], 'color','k','marker','.','markersize',16)
    set(gcf,'color','w');
    grid on
    xlabel("x [m]")
    ylabel("u [m/s]")

    subplot(3,1,2);
    plot([past_x x], [past_v v], 'color','k','marker','.','markersize',16)
    set(gcf,'color','w');
    grid on
    xlabel("x [m]")
    ylabel("v [m/s]")
    
    subplot(3,1,3);
    plot([past_x x], [past_w w], 'color','k','marker','.','markersize',16)
    set(gcf,'color','w');
    grid on
    xlabel("x [m]")
    ylabel("w [m/s]")
    
    %% 3D Trajectory
    % Convert color code to 1-by-3 RGB array (0~1 each)
    str = '#D9FFFF';
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

    str = '#fbd9d3';
    prev_color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    
    fig_num = fig_num + 1;
    figure(fig_num)
    h1 = plot3(x, y, z, '-o','Color','b','MarkerSize',10, 'MarkerFaceColor', color, 'DisplayName', 'Trajetória local otimizada');
    hold on
    h2 = plot3(past_x, past_y, past_z,'-o','Color','r','MarkerSize',10, 'MarkerFaceColor', prev_color, 'DisplayName', 'Trajetória já percorrida');
    h3 = plot3([x_d, past_x_d], [y_d, past_y_d], [z_d, past_z_d],'--','Color', [100, 100, 100]/255, 'DisplayName', 'Trajetória global');
    plot_obstacle(m, x_obs, y_obs, z_obs, R_obs)
    legend([h1 h2 h3])
    grid on
    axis equal
    set(gcf,'color','w');
    xlabel("x [m]")
    ylabel("y [m]")
    zlabel("z [m]")
    % view(-101.9, 1.2) % Different point of view
    hold off
    
    %% Performance and Control Limits 
    
    fig_num = fig_num + 1;
    figure(fig_num)
    plot(t_array(2:end), T, 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("T")
    set(gcf,'color','w');
    
    fig_num = fig_num + 1;
    figure(fig_num)
    plot(t_array, rad2deg(gamma), 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("\gamma")
    set(gcf,'color','w');
    
    fig_num = fig_num + 1;
    figure(fig_num)
    plot(t_array(2:end), rad2deg(phi), 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("\phi")
    set(gcf,'color','w');

    
    %% Historical Performance and Control Limits 

    fig_num = fig_num + 1;
    figure(fig_num)
    subplot(3,1,1);
    plot([past_x x(1:end-1)], [past_T T], 'color','k','marker','.','markersize',16)
    xlabel("x [m]")
    ylabel("T [N]")
    grid on
    
    subplot(3,1,2);   
    plot([past_x x], [rad2deg(past_gamma) rad2deg(gamma)], 'color','k','marker','.','markersize',16)
    xlabel("x [m]")
    ylabel("\gamma [º]")
    grid on
    
    subplot(3,1,3);
    plot([past_x x(1:end-1)], [rad2deg(past_phi) rad2deg(phi)], 'color','k','marker','.','markersize',16)
    xlabel("x [m]")
    ylabel("\phi [°]")
    grid on
    set(gcf,'color','w');
