function [e, S, b, T_max, mass, AR, k, C_D_0, ...
          gamma_max, gamma_min, phi_max, phi_min, ...
          n_max, n_min, V_min_drag, T_min] = vehicle_parameters(vehicle_case, G, RHO)

    switch vehicle_case
        case 0 % DECODE 2
            % Oswald factor [-]
            e = 0.8;

            % Wing area [m²]
            S = 1.554;

            % Span [m]
            b = 3.9422;

            % Maximum thrust (all operating engines) [N]
            T_max = 130;

            % Aircraft mass [kg]
            mass = 25;

            % Aspect Ratio [-]
            AR = b^2/S;

            % k value [-]
            k = 0.0334;

            % Cd0 value [-]
            % C_D_0 = 0.0715;
            C_D_0 = 0.0143;

            gamma_max = deg2rad(20);
            gamma_min = deg2rad(-20);

            phi_max = deg2rad(30);
            phi_min = deg2rad(-30);

            n_max = 3.8;
            n_min = -1.5;

            V_min_drag = sqrt(mass*G*2/(RHO*S))*((k/C_D_0)^(1/4));
            T_min = 1/2*RHO*S*2*C_D_0.*V_min_drag.^2;

        case 1 % DECODE 1
            % Oswald factor [-]
            e = 0.8;

            % Wing area [m²]
            S = 0.966;

            % Span [m]
            b = 2.9483;

            % Aircraft mass [kg]
            mass = 25;

            % Maximum thrust (all operating engines) [N]
            T_max = mass*G*0.361;

            % Aspect Ratio [-]
            AR = b^2/S;

            % k value [-]
            k = 0.0334;

            % Cd0 value [-]
            % C_D_0 = 0.0715;
            C_D_0 = 0.0143;

            gamma_max = deg2rad(20);
            gamma_min = deg2rad(-20);

            phi_max = deg2rad(40);
            phi_min = deg2rad(-40);

            n_max = 3.8;
            n_min = -1.5;

            V_min_drag = sqrt(mass*G*2/(RHO*S))*((k/C_D_0)^(1/4));
            T_min = 1/2*RHO*S*2*C_D_0.*V_min_drag.^2;
    end