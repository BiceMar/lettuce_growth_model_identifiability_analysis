modelName = 'lettuce_growth_model';

numDerivatives = 6;

syms c_alpha c_beta c_resp_sht c_gr_max c_epsilon_max 
parametersToAnalyze = [
    c_alpha; c_beta; c_resp_sht; c_gr_max
];

fprintf('Model: %s\n', modelName);
fprintf('Number of Lie derivatives: %d\n', numDerivatives);
disp('Parameters to analyze:');
disp(parametersToAnalyze.');

genssiMain(modelName, numDerivatives, parametersToAnalyze);
