# Identifiability of Lettuce Growth Model

This repository contains MATLAB scripts and results related to the analysis of a lettuce growth model. The analysis includes studying parameter correlations, sensitivity analysis, and structural identifiability of the model.

## Scripts description

* **`script/correlation.m`**: Analyzes the correlation between the parameters of the model. The results are visualized in `img/correlation.png`.
* **`script/sensitivity.m`**: Performs sensitivity analysis to understand how changes in model parameters affect the output. The results are visualized in `img/sensitivity.png`.
* **`script/structural_identifiability.m`**: Assesses the structural identifiability of the model using the Matlab GenSSI toolbox.


For a detailed discussion of the scripts and results, please refer to the document Result_discussion

## Dependencies

* MATLAB (Version 2024b)
* ([[GenSSI](https://github.com/genssi-developer/GenSSI)]) Toolbox
