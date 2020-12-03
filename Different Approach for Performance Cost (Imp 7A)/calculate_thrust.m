function T = calculate_thrust(u, v, w, ...
                              u_dot, v_dot, w_dot, ...
                              n, ...
                              mass, g, rho, S, C_D_0, k)

    %{
        Calculate thrust [N], given the speed [m/s] and 
        acceleration [m/s²] in each axis. The acceleration of 
        gravity [m/s²], the vehicle mass [kg], density of air [kg/m³],
        reference area [m²] and the dimensionless parameters 
        C_d_0 [-] and k [-] are also needed.
    %}
    V = sqrt(u^2 + v^2 + w^2);

    C_L = 2*n*mass*g/(rho*S*V^2);
    C_D = C_D_0 + k*C_L^2;

    D = 1/2*rho*S*C_D*V^2;

    T = (sqrt(w^2 + v^2 + u^2)*(mass*w*w_dot + g*mass*w + mass*v*v_dot + mass*u*u_dot) + D*(w^2 + v^2 + u^2))/(w^2 + v^2 + u^2);