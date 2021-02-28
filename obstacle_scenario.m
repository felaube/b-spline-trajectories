function [x_obs, y_obs, z_obs, R_obs, u_obs, v_obs, w_obs, num_obs] = obstacle_scenario(scenario_id)
    
    %{
        This function was written in order to store different scenario 
        for obstacles positioning. Given scenario_id, it returns
        corresponding initial positioning (x, y, z), radius (R),
        and velocities (u, v, w) of obstacles. In addition, it also returns
        the number of obstacles in the scenario.
    %}
    switch scenario_id
        case 0
            x_obs = [];
            y_obs = [];
            z_obs = [];
            R_obs = [0];
            u_obs = [];
            v_obs = [];
            w_obs = [];
        case 1 % Trajectory from Section V A of Alturbeh e Whidborne (2014)
            x_obs = 3500;
            y_obs = 1000;
            z_obs = 1000;
            R_obs = 200;
            u_obs = 0;
            v_obs = 0;
            w_obs = 0;
        case 1.5
            x_obs = 1000;
            y_obs = 1000;
            z_obs = 1000;
            R_obs = 200;
            u_obs = 0;
            v_obs = 0;
            w_obs = 0;
        case 2
            x_obs = 800;
            y_obs = 990;
            z_obs = 900;
            R_obs = 400;  
            u_obs = 0;
            v_obs = 0; 
            w_obs = 0;
        case 2.1
            x_obs = [800 1000];
            y_obs = [990 1690];
            z_obs = [900 900];
            R_obs = [400 150];
            u_obs = [0 0];
            v_obs = [0 0];
            w_obs = [0 0];
        case 2.2
            x_obs = [800 1000 3000];
            y_obs = [990 1690 1490];
            z_obs = [900 900 1000];
            R_obs = [400 150 200];
            u_obs = [0 0 0];
            v_obs = [0 0 0];
            w_obs = [0 0 0];
        case 3
            x_obs = 2000;
            y_obs = 1000;
            z_obs = 1000;
            R_obs = 200;
            u_obs = 0;
            v_obs = 0;
            w_obs = 0;
        case 4
            x_obs = [3500 4500];
            y_obs = [1000 990];
            z_obs = [1000 850];
            R_obs = [200 200];
            u_obs = [0 0];
            v_obs = [0 0];
            w_obs = [0 0];        
        case 5
            x_obs = 5000;
            y_obs = 1000;
            z_obs = 1000;
            R_obs = 200;
            u_obs = 0;
            v_obs = 0;
            w_obs = 0;
        case 6
            x_obs = 500;
            y_obs = 1000;
            z_obs = 1200;
            R_obs = 200;
            u_obs = 0;
            v_obs = 0;
            w_obs = 0;
        case 7
            x_obs = 6000;
            y_obs = 990;
            z_obs = 1000;
            R_obs = 400;
            u_obs = 0;
            v_obs = 0;
            w_obs = 0;
        case 8
            x_obs = [6000 8000 8750 10000];
            y_obs = [990 980 1000 990];
            z_obs = [1000 1200 900 1000];
            R_obs = [400 200 100 200];
            u_obs = [0 0 0 0];
            v_obs = [0 0 0 0];
            w_obs = [0 0 0 0];
        case 8.1 % Same as 8, but obstacles are closer to aircraft
            x_obs = [2000 4000 4750 6000];
            y_obs = [990 980 1000 990];
            z_obs = [1000 1200 900 1000];
            R_obs = [400 200 100 200];
            u_obs = [0 0 0 0];
            v_obs = [0 0 0 0];
            w_obs = [0 0 0 0];
        case 9 % Trajectory from Section V B of Alturbeh e Whidborne (2014)
            x_obs = [2000 2100];
            y_obs = [1000 1000];
            z_obs = [1000 1000];
            R_obs = [200 200];
            u_obs = [-18 15];
            v_obs = [0 0];
            w_obs = [0 0];
        case 10 % Like case 9, but the obstacle are static
            x_obs = [2000 2100];
            y_obs = [1000 1000];
            z_obs = [1000 1000];
            R_obs = [200 200];
            u_obs = [0 0];
            v_obs = [0 0];
            w_obs = [0 0];    
    end
    
    num_obs = length(x_obs);