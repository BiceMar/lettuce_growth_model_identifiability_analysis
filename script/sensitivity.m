clear; clc; close all;

fprintf('Sensitivity Analysis (Validity Checks Removed)\n');

% Load model parameters
params = load_parameters();

inputs.U_PAR = 150;  % W m^-2
inputs.U_CO2 = 600;  % ppm
inputs.U_T   = 20;   % degrees C

% Simulation time
t_start_days = 0;
t_end_days = 40;
t_span_seconds = [t_start_days, t_end_days] * 24 * 3600;

% Define evaluation time points 
num_time_points = 100;
t_eval_seconds = linspace(t_span_seconds(1), t_span_seconds(2), num_time_points);
t_eval_days = t_eval_seconds / (24 * 3600); % For plotting

% Initial conditions
X_nsdw_initial = 0.5; % g m^-2
X_sdw_initial  = 1.0; % g m^-2
X0 = [X_nsdw_initial; X_sdw_initial];

% ODE solver options
options = odeset('RelTol', 1e-6', 'AbsTol', 1e-8', 'NonNegative', [1, 2]);

% Nominal simulation
fprintf('Running nominal simulation...\n');
[~, X_nom] = ode45(@(t, X) lettuceODE(t, X, params, inputs), t_eval_seconds, X0, options);
DW_nom_ts = X_nom(:, 1) + X_nom(:, 2);
fprintf('Nominal Final DW = %.2f g m^-2\n', DW_nom_ts(end));

% Sensitivity analysis setup
parameters_to_analyze = fieldnames(params);
num_params = length(parameters_to_analyze);
perturb_factor = 0.05; % 5% perturbation

% Initialize matrices to store sensitivity time series
abs_sensitivity_matrix = NaN(num_time_points, num_params);
rel_sensitivity_matrix = NaN(num_time_points, num_params); % Matrix for relative sensitivity
scaled_sensitivity_matrix = NaN(num_time_points, num_params);

fprintf('\nStarting sensitivity analysis (perturbation = %.2f%%)...\n', perturb_factor * 100);

% Loop through parameters
for i = 1:num_params
    param_name = parameters_to_analyze{i};
    fprintf('Analyzing parameter: %s\n', param_name);

    nominal_value = params.(param_name);

    params_perturbed = params;
    delta_theta = nominal_value * perturb_factor;
    if delta_theta == 0
       fprintf('  Warning: delta_theta is zero for parameter %s.\n', param_name);
    end
    params_perturbed.(param_name) = nominal_value + delta_theta;

    [~, X_pert] = ode45(@(t, X) lettuceODE(t, X, params_perturbed, inputs), t_eval_seconds, X0, options);
    DW_pert_ts = X_pert(:, 1) + X_pert(:, 2);

    % Calculate sensitivity time series 

    % 1. Absolute sensitivity
    abs_sensitivity_ts = (DW_pert_ts - DW_nom_ts) / delta_theta;
    abs_sensitivity_matrix(:, i) = abs_sensitivity_ts;

    % 2. Relative sensitivity 
    rel_change_output = (DW_pert_ts - DW_nom_ts) ./ DW_nom_ts; 
    relative_sensitivity_ts = rel_change_output / perturb_factor;
    rel_sensitivity_matrix(:, i) = relative_sensitivity_ts;

    % 3. Scaled sensitivity
    mean_DW_nom = mean(DW_nom_ts);
    scale_factor_param = nominal_value;
    scale_factor_output = mean_DW_nom;
    if scale_factor_output == 0
         fprintf('  Warning: Mean nominal DW is zero, scaled sensitivity will be Inf/NaN.\n');
    end
    scaled_sensitivity_matrix(:, i) = abs_sensitivity_matrix(:, i) * (scale_factor_param / scale_factor_output);

end

fprintf('\nSensitivity analysis complete.\n');

% Visualization 
fprintf('\nPlotting time-dependent sensitivities...\n');

% Define line styles to cycle through
line_styles = {'-', '--', ':', '-.'};
num_styles = length(line_styles);

% Setup colors
colors_rel = lines(num_params);
legend_entries_rel = strrep(parameters_to_analyze, '_', '\_'); % Legend 

% Plot: Relative Sensitivity
figure;
hold on;
plot_handles_rel = [];
for k = 1:num_params 

    % colors
    current_color = colors_rel(mod(k-1, size(colors_rel, 1)) + 1, :); % Cycle colors

    % line styles
    current_style = line_styles{mod(k-1, num_styles) + 1};

    % Plot line
    h = plot(t_eval_days, rel_sensitivity_matrix(:, k), ...
             'LineWidth', 1.5, ...
             'Color', current_color, ...
             'LineStyle', current_style); 
    plot_handles_rel(k) = h;

end
hold off;
xlabel('Time (days)');
ylabel('Relative Sensitivity (rrSI_{k})');
title('Time-Dependent Relative Sensitivity of Total Dry Weigth');
legend(plot_handles_rel, legend_entries_rel, 'Location', 'best');

grid on;
ylim([-0.35 1.4]);
set(gcf, 'Position', [125, 125, 900, 600]);

% Calculate and plot importance indices 
fprintf('\nCalculating aggregate importance indices...\n');
rho_k_msqr = NaN(1, num_params);
rho_k_mabs = NaN(1, num_params);

mean_DW_nom_agg = mean(DW_nom_ts); 

for i = 1:num_params 
    try
        rho_k_msqr(i) = sqrt(mean(scaled_sensitivity_matrix(:, i).^2));
        rho_k_mabs(i) = mean(abs(scaled_sensitivity_matrix(:, i)));
    catch calc_err
        fprintf('Warning: Could not calculate rho_k for parameter %d: %s\n', i, calc_err.message);
        rho_k_msqr(i) = NaN;
        rho_k_mabs(i) = NaN;
    end
end

fprintf('Parameter Importance Indices (rho_k based on Scaled Sensitivity):\n');
fprintf('%-20s | %-10s | %-10s\n', 'Parameter', 'rho_msqr', 'rho_mabs');
fprintf('%-20s | %-10s | %-10s\n', repmat('-',1,20), repmat('-',1,10), repmat('-',1,10));

valid_rho_indices_agg = find(~isnan(rho_k_msqr));
if ~isempty(valid_rho_indices_agg)
    [sorted_rho_msqr, sort_idx_msqr] = sort(rho_k_msqr(valid_rho_indices_agg), 'descend'); % Sort valid ones
    sorted_names_agg = parameters_to_analyze(valid_rho_indices_agg(sort_idx_msqr));
    rho_mabs_sorted = rho_k_mabs(valid_rho_indices_agg(sort_idx_msqr)); % Get corresponding mabs

    for k_idx = 1:length(sort_idx_msqr)
         fprintf('%-20s | %-10.4f | %-10.4f\n', ...
             sorted_names_agg{k_idx}, ...
             sorted_rho_msqr(k_idx), ...
             rho_mabs_sorted(k_idx)); % Display sorted valid results
    end

    % Bar plot for rho_msqr (only plots the valid, sorted ones)
    figure;
    bar(sorted_rho_msqr);
    set(gca, 'XTick', 1:length(sorted_names_agg), 'XTickLabel', strrep(sorted_names_agg, '_', '\_'), 'XTickLabelRotation', 45);
    ylabel('Importance Index \rho_k^{msqr}');
    title('Parameter Importance Index (\rho_k^{msqr})');
    grid on;
    set(gcf, 'Position', [200, 200, 900, 600]);
else
    fprintf('No valid (non-NaN) rho_k indices were calculated to display/plot.\n');
end
