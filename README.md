# 2d_insilico_biopsy_ssa_validation

## Structure of this Repository

This is the code used in the puplication "2D reaction-diffusion model-based biopsy simulation for dynamic tumor growth parameter estimation" which is currently under peer review. It is sorted by the three main chapters of the paper's methodology (PDE solver, biopsy generation, validation), plus a script for the spectral-spatial analysis method in the version that has been used for the validation.

### PDE Solver

solution_RD_2D_exp and solution_RD_2D_log: MATLAB code to solve equations (1) and (2) with all tested initial conditions (Dirac-Delta approximations, Gaussian curve), on all tested grids. The initial condition and the grid can be selected inside the function definition.

The analytical solutions to the exponential growth equation are implemented as well.


### Biopsy Generation

biopsy_simulation: MATLAB code to generate the simulations, normalized exponential growth as demonstrated in the publication is pre-selected (outputs plots and optional excel-files with the coordinates of the cell positions).


### Spectral-Spatial Analysis

pattern_analysis: Python code with the method from Schlicke et al. (2026) to estimate diffusion and growth parameters from the biopsies (outputs npz-files with dictionaries of the estimates together with some additional information, as well as plots of the biopsies and plots of the data-2PCF and PSD together with the fitted theoretical curves). 

To execute this script, you need the script ssa_schlicke_25 which is currently not included in this repository. Please contact the authors for a personal copy.


### Validation

insilico_validation: Python code to perform the curve fits with the conversion laws.
