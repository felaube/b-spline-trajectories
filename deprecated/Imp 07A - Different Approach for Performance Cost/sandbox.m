t_array = linspace(0,200,51*200);
t_diff = diff(t_array);

u = 2*sin(t);
v = 2*cos(t/3).*t;

u_dot = 2*cos(t);
v_dot = 2*cos(t/3) - 2/3*sin(t/3).*t;

V = sqrt(u.^2 + v.^2);
V_dot_numerical = diff(V)./t_diff;
V_dot_analytical = sqrt(u_dot.^2 + v_dot.^2);

figure(1)
plot(t_array(2:end), V_dot_numerical, 'DisplayName', 'Numerical')
hold on
plot(t_array(2:end), V_dot_analytical(2:end), 'DisplayName', 'Analytical')

figure(2)
plot(t_array(2:end), (V_dot_analytical(2:end) - V_dot_numerical)./V_dot_analytical(2:end)*100)

