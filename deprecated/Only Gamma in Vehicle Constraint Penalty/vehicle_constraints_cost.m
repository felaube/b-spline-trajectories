function performance_margin = vehicle_constraints_cost(gamma)
   
% Aircraft: DECODE 2
    
gamma_max = deg2rad(20);
gamma_min = deg2rad(-20);

performance_margin = 1.5*(max(gamma)/gamma_max)^2;