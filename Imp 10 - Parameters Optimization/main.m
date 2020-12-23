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

opts2 = optimoptions('fmincon');
opts2.Display = 'off';
opts2.Algorithm = 'active-set';

% Scaling factors
% Default
LAMBDA_P = 1; % position
LAMBDA_S = 1; % speed
LAMBDA_PRF = 1; % vehicle contraints
LAMBDA_OB = 0; % obstacles constraints
LAMBDA_T = 1; % terminal attitude
LAMBDA_H = 1; % terminal heading angle
LAMBDA_F = 1; % terminal flight path angle

% Optimizin parameters
ic = [LAMBDA_P; LAMBDA_S; LAMBDA_PRF; LAMBDA_T];
objective = @(parameters) trajectory_optimization(parameters, opts);
lower_bounds = [1, 1, 1, 1];

tic
[optimalWayPoints, fval, ~, output] = fmincon(objective, ic(:), [],[],[],[],lower_bounds,[],[],opts2);
t = toc;

sprintf("Iterations: %i funcCount: %i", output.iterations, output.funcCount)
sprintf("Execution time: %.8f", t)
fprintf("For LAMBDA_P = %d, LAMBDA_S = %d, LAMBDA_PRF = %d and LAMBDA_T = %d, the smoothness was %f.\n", optimalWayPoints(1), optimalWayPoints(2), optimalWayPoints(3), optimalWayPoints(4), -fval)
