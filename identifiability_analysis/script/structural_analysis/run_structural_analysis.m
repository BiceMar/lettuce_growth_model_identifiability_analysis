
modelName = 'lettuce_growth_model';

% Set the maximum number of Lie derivatives to be computed by GenSSI.
% A higher number allows for deeper analysis but increases computation time.
numDerivatives = 6;

% Declare model parameters as symbolic variables.
% GenSSI performs symbolic computations, so parameters involved in the analysis,
% or present in the model equations, need to be defined symbolically.
syms c_alpha c_beta c_resp_sht c_gr_max c_epsilon_max

% Define a specific subset of parameters for structural identifiability analysis.
% Here, we focus on c_alpha, c_beta, c_resp_sht, and c_gr_max as explained in the report.
parametersToAnalyze = [
    c_alpha;      
    c_beta;       
    c_resp_sht;   
    c_gr_max      
];


fprintf('Model selected for analysis: %s\n', modelName);
fprintf('Maximum number of Lie derivatives to compute: %d\n', numDerivatives);
disp('Subset of parameters selected for identifiability analysis:');
disp(parametersToAnalyze.'); 

% Call the main GenSSI function to perform the structural identifiability analysis.
% This function will use the specified model, number of Lie derivatives, and
% the selected subset of parameters to determine their local and global
% structural identifiability properties.
genssiMain(modelName, numDerivatives, parametersToAnalyze);