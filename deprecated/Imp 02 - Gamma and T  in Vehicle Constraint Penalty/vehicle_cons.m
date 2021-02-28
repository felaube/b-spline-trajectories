function [c, ceq] = vehicle_cons(x, C_u, C_v, C_w, x_0, y_0, z_0, ... 
                                 R, R_dot, R_dot_dot, R_int, t_diff, ...
                                 x_d, y_d, z_d, psi_d, gamma_d, ...
                                 g, rho, e, S, b, T_max, mass, AR, k, C_D_0, T_min)
        
%% Construct vectors of coefficients 
C_u_ = [C_u(1:3); x(1:4)];
C_v_ = [C_v(1:3); x(5:8)];
C_w_ = [C_w(1:3); x(9:12)];

%% Calculate position
x = x_0 + C_u_'*R_int;
y = y_0 + C_v_'*R_int;
z = z_0 + C_w_'*R_int;
    
%% Calculate speeds and their corresponding derivatives
u = C_u_'*R;

v = C_v_'*R;

w = C_w_'*R;

V = sqrt(u.^2 + v.^2 + w.^2);

%% Calculate flight path angle, heading angle and their derivatives
gamma = asin(w./V);
psi = asin(v./(V.*cos(gamma)));

gamma(isnan(gamma)) = 0; % Not sure if this is necessary. CHECK LATER
psi(isnan(psi)) = 0; % Not sure if this is necessary. CHECK LATER

psi_dot = diff(psi)./t_diff;

gamma_dot = diff(gamma)./t_diff;
gamma_dot_dot = diff(gamma_dot)./t_diff(2:end);

%% Calculate bank angle
phi = atan2(psi_dot(1:end).*V(2:end).*cos(gamma(2:end)), g*cos(gamma(2:end)) + V(2:end).*gamma_dot);

c = 0;

ceq(1) = gamma(end) - gamma_d(end);
ceq(2) = psi(end) - psi_d(end);
ceq(3) = phi(end);
ceq(4) = x(end) - x_d(end);
ceq(5) = y(end) - y_d(end);
ceq(6) = z(end) - z_d(end);
