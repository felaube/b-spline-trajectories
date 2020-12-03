function phi = calculate_phi(u, v, w, ...
                             u_dot, v_dot, w_dot, ...
                             u_dot_dot, v_dot_dot, w_dot_dot, g)

    %{
        Calculate bank angle [rad], given the speed [m/s] and 
        acceleration [m/s²] in each axis. The acceleration of 
        gravity [m/s²] is also needed.
    %}


    phi = atan2(((u^2*(v_dot)-u*(u_dot)*v)*sqrt(w^2+v^2+u^2)), ((abs(u)*v^2+u^2*abs(u))*(w_dot)+(-abs(u)*v*(v_dot)-u*abs(u)*(u_dot))*w+g*abs(u)*v^2+g*u^2*abs(u)));