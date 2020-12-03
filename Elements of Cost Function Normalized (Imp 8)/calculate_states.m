function [u, v, w, u_dot, v_dot, w_dot] = calculate_states(phi, T, n, ...
                                                           gamma, psi, ...
                                                           V, D, ...
                                                           mass, g)
    %{
        Calculate speed and acceleration in each axis, given the current
        bank angle [rad], thrust [N], load factor [-], flight path angle
        (which in this case is equal to the angle of attack) [rad] and
        the heading angle [rad]. The total speed [m/s], drag [N], vehicle 
        mass [kg] and acceleration of gravity [m/s²] are also needed.
    %}
    u = V*cos(gamma)*cos(psi);
    v = V*cos(gamma)*sin(psi);
    w = V*sin(gamma);
    
    V_dot = (T - D)/mass - g*sin(gamma);
    gamma_dot = (g/V)*(n*cos(phi) - cos(psi));
    psi_dot = (g/V)*(n*sin(phi)/cos(gamma));
    
    w_dot = (gamma_dot*V^2*(1-(w/V)^2)+w*V_dot)/V;
    v_dot = (psi_dot*(1-(v/(V*cos(gamma)))^2)*V^2*(cos(gamma))^2+v*(V_dot*cos(gamma)-V*sin(gamma)*gamma_dot))/(V*cos(gamma));
    u_dot = (V_dot*V - w_dot*w - v_dot*v)/u;

