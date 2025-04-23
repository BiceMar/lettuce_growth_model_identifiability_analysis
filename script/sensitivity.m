clear; clc; close all;

fprintf('Setting up Lettuce Growth Model Sensitivity Analysis...\n');

% Load model parameters
params = load_parameters();

inputs.U_PAR = 150;  % W m^-2 
inputs.U_CO2 = 600;  % ppm
inputs.U_T   = 20;   % degrees C

% Simulation time
t_start_days = 0;
t_end_days = 40; 
t_span_seconds = [t_start_days, t_end_days] * 24 * 3600; 

% Initial conditions
X_nsdw_initial = 0.5; % g m^-2
X_sdw_initial  = 1.0; % g m^-2
X0 = [X_nsdw_initial; X_sdw_initial];

% ODE solver options
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NonNegative', [1, 2]); 

% Simulate the system from start to end time
%  t_nom-> column vector of time points that the solver used internally
%  X_nom-> matrix of the solution values at the time points in t_nom.
fprintf('Running nominal simulation with updated parameters...\n');
[t_nom, X_nom] = ode45(@(t, X) lettuceODE(t, X, params, inputs), t_span_seconds, X0, options);

% Final total dry weight
X_nsdw_nom_final = X_nom(end, 1);
X_sdw_nom_final = X_nom(end, 2);
DW_nom_final = X_nsdw_nom_final + X_sdw_nom_final;
fprintf('Nominal Final DW = %.2f g m^-2\n', DW_nom_final);

% Sensitivity analysis setup
parameters_to_analyze = fieldnames(params); 
perturb_factor = 0.10; % +/- 10% perturbation
sensitivity_results = struct();

fprintf('\nStarting Sensitivity Analysis (Perturbation = %.1f%%)...\n', perturb_factor * 100);

% Loop through parameters and perfrorm forward difference approximation for the sensitivity
for i = 1:length(parameters_to_analyze)
    param_name = parameters_to_analyze{i};
    fprintf('Analyzing parameter: %s\n', param_name);

    nominal_value = params.(param_name);

    % Positive perturbation
    params_perturbed_pos = params;
    params_perturbed_pos.(param_name) = nominal_value * (1 + perturb_factor);

    % simulation with positive perturbation
    try
        [~, X_pert_pos] = ode45(@(t, X) lettuceODE(t, X, params_perturbed_pos, inputs), t_span_seconds, X0, options);
        DW_pert_pos_final = X_pert_pos(end, 1) + X_pert_pos(end, 2);
    catch ME
        warning('Simulation failed for + perturbation of %s: %s', param_name, ME.message);
        DW_pert_pos_final = NaN; 
    end

    % Calculate sensitivity
    if ~isnan(DW_pert_pos_final) && DW_nom_final ~= 0
        relative_change_output_pos = (DW_pert_pos_final - DW_nom_final) / DW_nom_final;
        relative_change_input_pos = perturb_factor;
        sensitivity_pos = relative_change_output_pos / relative_change_input_pos;
    elseif DW_nom_final == 0
         sensitivity_pos = NaN; 
    else
        sensitivity_pos = NaN;
    end

    sensitivity_results.(param_name) = sensitivity_pos;

end

fprintf('\nSensitivity Analysis Complete.\n');
fprintf('\nSensitivity Results (Relative Change in Final DW / Relative Change in Parameter):\n');

param_names = fieldnames(sensitivity_results); 
num_params = length(param_names);
sens_values = NaN(1, num_params); % allocate vector with NaNs

for k = 1:num_params
    value = sensitivity_results.(param_names{k});
    sens_values(k) = value; % Assign to the vector
end

% Process valid results
% Remove NaN values
valid_indices = ~isnan(sens_values);
param_names_valid = param_names(valid_indices); % Filter names
sens_values_valid = sens_values(valid_indices); % Filter values

% Sort by absolute sensitivity
if ~isempty(sens_values_valid) 
    [~, sort_idx] = sort(abs(sens_values_valid), 'descend'); 
    sorted_sens_values = sens_values_valid(sort_idx);
    sorted_param_names = param_names_valid(sort_idx); 

    % Display sorted results
    for i = 1:length(sorted_sens_values)
        fprintf('%-15s: %+.4f\n', sorted_param_names{i}, sorted_sens_values(i));
    end

    % Bar plot
    figure;
    bar(sorted_sens_values);
    set(gca, 'XTick', 1:length(sorted_param_names), 'XTickLabel', sorted_param_names, 'XTickLabelRotation', 45);
    ylabel('Relative Sensitivity');
    title(['Sensitivity of Final Total DW to Model Parameters (', num2str(perturb_factor*100), '% Perturbation)']);
    grid on;
    set(gcf, 'Position', [100, 100, 900, 600]); 
    set(gca, 'TickLabelInterpreter', 'none'); 

else
    fprintf('No valid sensitivity values found to display or sort.\n');
   
end