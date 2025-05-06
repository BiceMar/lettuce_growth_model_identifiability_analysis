clear; clc; close all;

fprintf('Lettuce Growth Model Correlation analysis\n');

% Load model parameters
params = load_parameters();

parameters_to_analyze = fieldnames(params); 
num_params = length(parameters_to_analyze);

% Environmental inputs
inputs.U_PAR = 150;  % W m^-2 
inputs.U_CO2 = 600;  % ppm
inputs.U_T   = 20;   % degrees C

% Simulation time
t_start_days = 0;
t_end_days = 40; 
t_span_seconds = [t_start_days, t_end_days] * 24 * 3600; 

% Time points for sensitivity evaluation
num_time_points = 100; % Number of points to evaluate sensitivity at
t_eval = linspace(t_span_seconds(1), t_span_seconds(end), num_time_points)'; 

% Initial conditions
X_nsdw_initial = 0.5; % g m^-2
X_sdw_initial  = 1.0; % g m^-2
X0 = [X_nsdw_initial; X_sdw_initial];

% ODE solver options
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NonNegative', [1, 2]);

% Run simulation at t_eval points
fprintf('Running nominal simulation for sensitivity analysis...\n');

[~, X_nom_traj] = ode45(@(t, X) lettuceODE(t, X, params, inputs), t_eval, X0, options);
DW_nom_traj = X_nom_traj(:,1) + X_nom_traj(:,2); % Sun ot dw and nsdw = total dw 

% Correlation Method
fprintf('\nStarting Correlation Method...\n');

% Calculate sensitivity matrix S using finite differences
S = zeros(num_time_points, num_params); 
step_fraction = 1e-4; % Relative step size for finite differences
fprintf('Calculating Sensitivity Matrix S (using Central Differences)...\n');
for j = 1:num_params
    param_name = parameters_to_analyze{j};
    fprintf('  Calculating sensitivity for: %s\n', param_name);

    nominal_value = params.(param_name);
    delta_param = nominal_value * step_fraction;

    if delta_param ~= 0
        % Positive perturbation
        params_perturbed_pos = params;
        params_perturbed_pos.(param_name) = nominal_value + delta_param;

        try
            [~, X_pert_traj_pos] = ode45(@(t, X) lettuceODE(t, X, params_perturbed_pos, inputs), t_eval, X0, options);
            DW_pert_traj_pos = X_pert_traj_pos(:,1) + X_pert_traj_pos(:,2);
        catch ME
            warning('Simulation failed for + perturbation of %s: %s', param_name, ME.message);
            DW_pert_traj_pos = NaN(num_time_points, 1);
        end

        % Negative perturbation
        params_perturbed_neg = params;
        params_perturbed_neg.(param_name) = nominal_value - delta_param;

        try
            [~, X_pert_traj_neg] = ode45(@(t, X) lettuceODE(t, X, params_perturbed_neg, inputs), t_eval, X0, options);
            DW_pert_traj_neg = X_pert_traj_neg(:,1) + X_pert_traj_neg(:,2);
        catch ME
            warning('Simulation failed for - perturbation of %s: %s', param_name, ME.message);
            DW_pert_traj_neg = NaN(num_time_points, 1);
        end

        % Comutation of central difference sensitivity
        if all(~isnan(DW_pert_traj_pos)) && all(~isnan(DW_pert_traj_neg))
            S(:, j) = (DW_pert_traj_pos - DW_pert_traj_neg) / (2 * delta_param);
        else
            S(:, j) = NaN;
        end
    else
        warning("warning: delta_param = 0");
        S(:, j) = NaN;
    end
end

fprintf('Sensitivity matrix S calculated.\n');

% Remove columns with NaN (failed simulations)
valid_sens_cols = ~all(isnan(S), 1) & ~all(S == 0, 1);
S_valid = S(:, valid_sens_cols);
valid_param_names = parameters_to_analyze(valid_sens_cols);
num_valid_params = size(S_valid, 2);

if isempty(S_valid)
    error('Sensitivity matrix calculation failed for all parameters.');
end

% Rank check
fprintf('\n--- Rank Analysis ---\n');
StS = S_valid' * S_valid;
rank_StS = rank(StS); 

fprintf('Number of parametersof the model: %d\n', num_params);
fprintf('Number of parameters considered: %d\n', num_valid_params);
fprintf('Rank of S^T * S: %d\n', rank_StS);

if rank_StS < num_valid_params
    fprintf('Result: S^T*S is RANK DEFICIENT. Parameters are locally non-identifiable.\n');
else
    fprintf('Result: S^T*S is FULL RANK. Parameters are locally identifiable.\n');
end

% Pairwise correlation analysis
fprintf('\n--- Pairwise Correlation Analysis ---\n');
if num_valid_params > 1
    R_sensitivity = corrcoef(S_valid); % Correlation between columns of matrix S_valid

    % Display correlation matrix heatmap
    figure;
    imagesc(abs(R_sensitivity)); % Show absolute correlations
    colorbar;
    title('Absolute correlation between sensitivity columns');
    set(gca, 'XTick', 1:num_valid_params, 'XTickLabel', valid_param_names, 'XTickLabelRotation', 45);
    set(gca, 'YTick', 1:num_valid_params, 'YTickLabel', valid_param_names);
    set(gca, 'TickLabelInterpreter', 'none'); 
    axis square;

    % Report highly correlated pairs
    %correlation_threshold = 0.95; % threshold
    %fprintf('Pairs with absolute correlation >= %.2f:\n', correlation_threshold);
    %found_high_corr = false;
    %for i = 1:num_valid_params
    %    for j = i+1:num_valid_params 
    %        if abs(R_sensitivity(i, j)) >= correlation_threshold
    %            fprintf('  %s <--> %s : %.4f\n', valid_param_names{i}, valid_param_names{j}, R_sensitivity(i, j));
    %            found_high_corr = true;
    %        end
    %    end
    %end
    %if ~found_high_corr
    %    fprintf('  None found.\n');
    %end
else
    fprintf('Only one valid parameter sensitivity column. Cannot calculate correlations.\n');
end

% Total Correlation Analysis
fprintf('\n--- Total Correlation Analysis ---\n');
if num_valid_params > 1
    delta_tc = 0.95; % Threshold for total correlation calculation
    TC = zeros(1, num_valid_params);
    abs_R = abs(R_sensitivity);
    for i = 1:num_valid_params
        for j = 1:num_valid_params
            if i ~= j && abs_R(i, j) >= delta_tc
                TC(i) = TC(i) + abs_R(i,j)^2; % Sum of squared high correlations
                
            end
        end
    end

    % Sort by total correlation
    [sorted_TC, sort_idx_tc] = sort(TC, 'descend');
    sorted_TC_params = valid_param_names(sort_idx_tc);

    fprintf('Total Correlation (Sum of R^2 for |R|>=%.2f):\n', delta_tc);
    for i = 1:num_valid_params
        fprintf('  %-15s: %.4f\n', sorted_TC_params{i}, sorted_TC(i));
    end

    if sorted_TC(1) > 0
       fprintf('Parameter with highest total correlation: %s (TC = %.4f)\n', sorted_TC_params{1}, sorted_TC(1));
    else
       fprintf('No parameters found with high correlations above threshold.\n');
    end
else
     fprintf('Cannot calculate total correlations with less than 2 parameters.\n');
end

fprintf('\n Analysis Complete.\n');

