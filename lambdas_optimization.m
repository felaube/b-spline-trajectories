%% Prepare workspace
close all
clear

%% User Defined Parameters
% Results spreadsheet
SPREADSHEET = 'response.xls';
writematrix(["LAMBDA_P" "LAMBDA_S" "LAMBDA_LIM" "LAMBDA_AT" "convergence" "smoothness"], SPREADSHEET)

% Optimization options

opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

%% Hypercube
v = 1:10;

PARAMETERS = nchoosek(v,4);
for row = 1 : size(PARAMETERS, 1) 
    permutated_rows = perms(PARAMETERS(row, :));
    for perm_row = 1 : size(permutated_rows, 1)
        LAMBDA_P = permutated_rows(perm_row, 1); % position
        LAMBDA_S = permutated_rows(perm_row, 2); % speed
        LAMBDA_LIM = permutated_rows(perm_row, 3); % vehicle contraints
        LAMBDA_AT = permutated_rows(perm_row, 4); % terminal attitude
        
        lambdas = [LAMBDA_P LAMBDA_S LAMBDA_LIM LAMBDA_AT];
        trajectory_evaluation = trajectory_optimization(lambdas, opts);
        
        convergence = trajectory_evaluation(1); 
        smoothness = trajectory_evaluation(2);
        
        writematrix([LAMBDA_P LAMBDA_S LAMBDA_LIM LAMBDA_AT convergence smoothness], SPREADSHEET, 'WriteMode', 'append')
    end
end