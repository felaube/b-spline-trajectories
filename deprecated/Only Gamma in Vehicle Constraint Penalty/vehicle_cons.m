function [c, ceq] = vehicle_cons(x, C_u, C_v, C_w, R, R_dot, R_dot_dot, R_int, t_diff, x_d, y_d, z_d, psi_d, gamma_d)

ceq = 0;

%% Construct vectors of coefficients 
C_u_ = [C_u(1:3); x(1:4)];
C_v_ = [C_v(1:3); x(5:8)];
C_w_ = [C_w(1:3); x(9:12)];

%% Calculate position
x = 0 + C_u_'*R_int;
y = 10 + C_v_'*R_int;
z = 1200 + C_w_'*R_int;
    
%% Calculate speeds and their corresponding derivatives
u = C_u_'*R;
% u_dot = C_u_'*R_dot*1/t_h;
% u_dot_dot = C_u_'*R_dot_dot*(1/(t_h^2));

v = C_v_'*R;
% v_dot = C_v_'*R_dot*1/t_h;
% v_dot_dot = C_v_'*R_dot_dot*(1/(t_h^2));

w = C_w_'*R;
% w_dot = C_w_'*R_dot*1/t_h;
% w_dot_dot = C_w_'*R_dot_dot*(1/(t_h^2));

V = sqrt(u.^2 + v.^2 + w.^2);
V_dot = diff(V)./t_diff;

%% Calculate flight path angle, heading angle and their derivatives
gamma = asin(w./V);
psi = asin(v./(V.*cos(gamma)));

gamma(isnan(gamma)) = 0; % Not sure if this is necessary. CHECK LATER
psi(isnan(psi)) = 0; % Not sure if this is necessary. CHECK LATER

psi_dot = diff(psi)./t_diff;
% psi_dot_dot = diff(psi_dot)./t_diff(2:end);

gamma_dot = diff(gamma)./t_diff;
gamma_dot_dot = diff(gamma_dot)./t_diff(2:end);
    
gamma_max = deg2rad(20);
gamma_min = deg2rad(-20);

%c(1:length(gamma)) = gamma - gamma_max;
%c(1 + length(gamma):2*length(gamma)) = gamma_min - gamma;

 % Aircraft: DECODE 2

%% Constants 
% Oswald factor [-]
% e = 0.8;

% Wing area [m²]
S = 1.554;

% Span [m]
b = 3.9422;

% Maximum thrust (all operating engines) [N]
T_max = 130;

% Aircraft mass [kg]
m = 25;

% Gravitational acceleration [m/s²]
g = 9.81;

% Air density (speed is considered to be AES) [kg/m³]
rho = 1.225;

% Aspect Ratio [-]
%A = b^2/S;

% k value [-]
k = 0.0334;

% Cd0 value (guessed) [-]
C_D_0 = 0.0715;

%% Calculate bank angle
phi = atan2(psi_dot(2:end).*V(3:end).*cos(gamma(3:end)), g*cos(gamma(3:end)) + V(3:end).*gamma_dot_dot);

%% Calculate load factor
n = (g*cos(gamma(3:end)) + V(3:end).*gamma_dot(2:end))./(g*cos(phi));

%% Calculate thrust
C_L = 2*n*m*g./(rho*S*V(3:end).^2);
C_D = C_D_0 + k*C_L.^2;

D = 1/2*rho*S*C_D.*V(3:end).^2;

T = D + m*V_dot(2:end) + m*g*sin(gamma(3:end));

T(isnan(T)) = 0; % Not sure if this is necessary. CHECK LATER

V_min_drag = sqrt(m*g*2/(rho*S))*((k/C_D_0)^(1/4));
T_min = 1/2*rho*S*2*C_D_0.*V_min_drag.^2;

c(1 : length(T)) = T - T_max;
c(length(T) + 1 : 2*length(T)) = T_min - T;

%c(2*length(gamma) + length(T) + 1) = x_d(end) - x(end);
%c(2*length(gamma) + length(T) + 2) = y_d(end) - y(end);
%c(2*length(gamma) + length(T) + 3) = z_d(end) - z(end);
ceq(1) = gamma(end) - gamma_d(end);
ceq(2) = psi(end) - psi_d(end);
ceq(3) = phi(end);
ceq(4) = x(end) - x_d(end);
ceq(5) = y(end) - y_d(end);
ceq(6) = z(end) - z_d(end);

% disp(0);
%c(:) = 0;
%ceq(:) = 0;