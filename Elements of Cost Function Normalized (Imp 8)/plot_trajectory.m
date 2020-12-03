function [x, y, z, u, v, w, ...
          u_dot, v_dot, w_dot, ...
          u_dot_dot, v_dot_dot, w_dot_dot, ...
          V, V_dot, ...
          gamma, psi, psi_dot, gamma_dot, ...
          psi_dot_dot, gamma_dot_dot, ...
          phi, n, C_L, C_D, D, T] = plot_trajectory(x_0, y_0, z_0, C_u, C_v, C_w, ...
                                                    R_int, R, R_dot, R_dot_dot, ...
                                                    t_array, t_diff, t_h, ...
                                                    g, rho, S, C_D_0, k, mass, ...
                                                    x_obs, y_obs, z_obs, R_obs, m, ...
                                                    past_x, past_y, past_z, ...
                                                    x_d, y_d, z_d)

    x = x_0 + C_u'*R_int;
    y = y_0 + C_v'*R_int;
    z = z_0 + C_w'*R_int;
    u = C_u'*R;
    v = C_v'*R;
    w = C_w'*R;
    u_dot = C_u'*R_dot/t_h;
    v_dot = C_v'*R_dot/t_h;
    w_dot = C_w'*R_dot/t_h;
    u_dot_dot = C_u'*R_dot_dot/(t_h^2);
    v_dot_dot = C_v'*R_dot_dot/(t_h^2);
    w_dot_dot = C_w'*R_dot_dot/(t_h^2);

    V = sqrt(u.^2 + v.^2 + w.^2);
    V_dot = diff(V)./t_diff;

    gamma = asin(w./V);
    psi = asin(v./(V.*cos(gamma)));

    psi_dot = diff(psi)./t_diff;
    psi_dot_dot = diff(psi_dot)./t_diff(2:end);

    gamma_dot = diff(gamma)./t_diff;
    gamma_dot_dot = diff(gamma_dot)./t_diff(2:end);


    %% Calculate bank angle
    phi = atan2(psi_dot(1:end).*V(2:end).*cos(gamma(2:end)), g*cos(gamma(2:end)) + V(2:end).*gamma_dot);

    %% Calculate load factor
    n = (g*cos(gamma(2:end)) + V(2:end).*gamma_dot(1:end))./(g*cos(phi));

    %% Calculate thrust
    C_L = 2*n*mass*g./(rho*S*V(2:end).^2);
    C_D = C_D_0 + k*C_L.^2;

    D = 1/2*rho*S*C_D.*V(2:end).^2;

    T = D + mass*V_dot(1:end) + mass*g*sin(gamma(2:end));
    
    %{
    figure(1)
    subplot(2,1,1);
    plot(x, y, 'color','k','marker','.','markersize',16)
    for i = 1 : m
        viscircles([x_obs(i) y_obs(i)], R_obs(i));
    end
    grid on
    axis equal
    xlabel("x")
    ylabel("y")

    %figure(2)
    subplot(2,1,2);
    plot(x, z, 'color','k','marker','.','markersize',16)
    for i = 1 : m
        viscircles([x_obs(i) z_obs(i)], R_obs(i));
    end
    grid on
    axis equal
    xlabel("x")
    ylabel("z")

    %{
    figure(3)
    plot(y, z, 'color','k','marker','.','markersize',16)
    xlabel("y")
    ylabel("z")
    %}

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
    % Convert color code to 1-by-3 RGB array (0~1 each)
    str = '#D9FFFF';
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    
    str = '#FBD9D3';
    prev_color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    
    figure(8)
    plot3(x, y, z,'-o','Color','b','MarkerSize',10, 'MarkerFaceColor', color)
    hold on
    plot3(past_x, past_y, past_z,'-o','Color','r','MarkerSize',10, 'MarkerFaceColor', prev_color)
    plot3(x_d, y_d, z_d,'-o','Color',[75 75 75]/255, 'MarkerSize',10, 'MarkerFaceColor', [200 200 200]/255)
    plot_obstacle(m, x_obs, y_obs, z_obs, R_obs)
    grid on
    axis equal
    xlabel("x")
    ylabel("y")
    zlabel("z")
    %view(-101.9, 1.2)
    hold off

    %{
    figure(9)
    plot(t_array, z, 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("z")
    %}
    
    %{
    figure(10)
    %plot(t_array, [1 1 n], 'color','k','marker','.','markersize',16)
    plot(t_array, [1 n], 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("n")
    %}
    
    %%{
    figure(11)
    %plot(t_array(3:end), T, 'color','k','marker','.','markersize',16)
    plot(t_array(2:end), T, 'color','k','marker','.','markersize',16)
    xlabel("t")
    ylabel("T")
    %%}