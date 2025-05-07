# Identifiability and Parameter Estimation of a Lettuce Growth Model

This repository contains analysis related to the lettuce growth model described in *"Validation of a dynamic lettuce growth model for greenhouse climate control"* by Van Henten, in particular the identifiability analysis and the parameter estimation.

1. **Identifiability Analysis:**
    * **`identifiability_analysis/sensitivity_analysis_scripts/sensitivity.m`**: runs sensitivity analysis. 
    * **`identifiability_analysis/sensitivity_analysis_scripts/correlation.m`**: runs the correlation method. 
    * **`identifiability_analysis/structural_analysis_scripts/run_structural_analysis.m`**: runs the structural identifiability analysis of the model using the GenSSI toolbox.
    
    Results are visualized in `identifiability_analysis/results`. For a detailed discussion of these MATLAB scripts and results, please refer to the document `identifiability_analysis/report_identifiability_analysis.pdf`.

2. **Parameter Estimation:**
   *  **`parameter_estimation/param_estimation_simulated_data.ipynb`**: Jupyter Notebook focused on estimating model parameters from simulated noisy data, including post-estimation analysis such as standard errors and metrics to measure the goodness of the fit.


## Versions and Tools

* **MATLAB:** Version 2024b 
    * [GenSSI Toolbox](https://github.com/genssi-developer/GenSSI)
* **Python:** Version 3.11.0
