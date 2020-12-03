function n = calculate_n(u, v, w, ...
                         u_dot, v_dot, w_dot, ...
                         u_dot_dot, v_dot_dot, w_dot_dot, ...
                         g)
    %{ 
        Calculate load factor, given the speed [m/s] and acceleration [m/s²] 
        in each axis. The acceleration of gravity [m/s²] is also needed.
    %}

    n = (((v^2+u^2)*(w_dot)+(-v*(v_dot)-u*(u_dot))*w+g*v^2+g*u^2)*sqrt(((v^4+2*u^2*v^2+u^4)*(w_dot)+(((-2*v^3-2*u^2*v)*(v_dot)-2*u*(u_dot)*v^2-2*u^3*(u_dot))*w+2*g*v^4 + 4*g*u^2*v^2+2*g*u^4)*(w_dot)+((v^2+u^2)*(v_dot)^2+(u_dot)^2*v^2+u^2*(u_dot)^2)*w^2+((-2*g*v^3-2*g*u^2*v)*(v_dot)-2*g*u*(u_dot)*v^2-2*g*u^3*(u_dot))*w+(u^2*v^2+u^4)*(v_dot)^2+(-2*u*(u_dot)*v^3-2*u^3*(u_dot)*v)*(v_dot)+((u_dot)^2+g^2)*v^4+(u^2*(u_dot)^2+2*g^2*u^2)*v^2+g^2*u^4)/((v^4+2*u^2*v^2+u^4)*(w_dot)^2+(((-2*v^3-2*u^2*v)*(v_dot)-2*u*(u_dot)*v^2-2*u^3*(u_dot))*w+2*g*v^4+4*g*u^2*v^2+2*g*u^4)*(w_dot)+(v^2*(v_dot)^2+2*u*(u_dot)*v*(v_dot)+u^2*(u_dot)^2)*w^2+((-2*g*v^3-2*g*u^2*v)*(v_dot)-2*g*u*(u_dot)*v^2-2*g*u^3*(u_dot))*w+g^2*v^4+2*g^2*u^2*v^2+g^2*u^4)))/(g*sqrt(v^2+u^2)*sqrt(w^2+v^2+u^2));