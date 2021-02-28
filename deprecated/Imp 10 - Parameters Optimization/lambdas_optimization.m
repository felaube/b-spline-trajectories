%% Prepare workspace
close all
clear

%% User Defined Parameters

% Optimization options

%opts = optimoptions('fmincon', 'UseParallel', true);
opts = optimoptions('fmincon');
%opts.Display = 'iter';
opts.Display = 'off';
opts.Algorithm = 'sqp';
%opts.Algorithm = 'interior-point';
%opts.Algorithm = 'active-set';
%opts.MaxFunEvals = 5000;
%opts.MaxFunEvals = 10000;
%opts.MaxFunEvals = inf;
%opts.MaxIterations = inf;
%opts.StepTolerance = 1e-16;
%opts.TolX = 1e-16;
%opts.ConstraintTolerance = 1e-03;

%opts2 = optimoptions('fmincon');
%opts2.Display = 'off';
%opts2.Algorithm = 'sqp';

% opts2 = optimset('Display', 'off', 'TolX', 0.1);
opts2 = optimoptions('fmincon', ...
                     'Algorithm', 'sqp', ...
                     'Display', 'off', ...
                     'TolX', 0.1, ...
                     'TolFun', 0.01);

% Scaling factors
% Default
LAMBDA_P = 1; % position
LAMBDA_S = 1; % speed
LAMBDA_PRF = 1; % vehicle contraints
LAMBDA_OB = 0; % obstacles constraints
LAMBDA_T = 1; % terminal attitude
LAMBDA_H = 1; % terminal heading angle
LAMBDA_F = 1; % terminal flight path angle

%% Optimizing parameters
%{
ic = [LAMBDA_P; LAMBDA_S; LAMBDA_PRF; LAMBDA_T];
objective = @(parameters) trajectory_optimization(parameters, opts);

tic
[optimalWayPoints, fval, ~, output] = fmincon(objective, ic(:), [],[],[],[],[],[],[],opts2);
% [optimalWayPoints, fval, ~, output] = fminsearch(objective, ic(:), opts2);
t = toc;

sprintf("Iterations: %i funcCount: %i", output.iterations, output.funcCount)
sprintf("Execution time: %.8f", t)
fprintf("For LAMBDA_P = %d, LAMBDA_S = %d, LAMBDA_PRF = %d and LAMBDA_T = %d. Objective: %f.\n", optimalWayPoints(1), optimalWayPoints(2), optimalWayPoints(3), optimalWayPoints(4), fval)
%}

%% Latin Hypercube
v = 1:10;

PARAMETERS = nchoosek(v,4);
for row = 1 : size(PARAMETERS, 1) 
    permutated_rows = perms(PARAMETERS(row, :));
    for perm_row = 1 : size(permutated_rows, 1)
        LAMBDA_P = permutated_rows(perm_row, 1); % position
        LAMBDA_S = permutated_rows(perm_row, 2); % speed
        LAMBDA_PRF = permutated_rows(perm_row, 3); % vehicle contraints
        LAMBDA_T = permutated_rows(perm_row, 4); % terminal attitude
        
        parameters = [LAMBDA_P LAMBDA_S LAMBDA_PRF LAMBDA_T];
        trajectory_evaluation = trajectory_optimization(parameters, opts);
        
        energy = trajectory_evaluation(1);
        convergence = trajectory_evaluation(2); 
        smoothness = trajectory_evaluation(3);
        negative_T = trajectory_evaluation(4);
         
        if (row == 1 && perm_row == 1)
            writematrix(["LAMBDA_P" "LAMBDA_S" "LAMBDA_PRF" "LAMBDA_T" "energy" "convergence" "smoothness" "negative_T"], 'response.xls')
        end
        
        writematrix([LAMBDA_P LAMBDA_S LAMBDA_PRF LAMBDA_T energy convergence smoothness, negative_T], 'response.xls', 'WriteMode', 'append')
    end
end