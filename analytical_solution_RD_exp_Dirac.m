% Analytical solution of the dim-dimensional reaction-diffusion equation with 
% exponential growth and Dirac-Delta initial condition.

function [u, r] = analytical_solution_RD_exp_Dirac(dim, t, rmax, gamma, D, delta_r)
    % grid
    Nr = rmax / delta_r;
    r = linspace(0, rmax, Nr+1)';

    % solution
    u = 1 / (4 * pi * D * t)^(dim/2) .* exp(-r.^2/(4 * D * t) + gamma * t);
end