% Analytical solution of the dim-dimensional reaction-diffusion equation with 
% exponential growth and a Gaussian curve as initial condition.

function u = analytical_solution_RD_exp_Gauss(dim, t, rmax, gamma, D, delta_r)
    % delta_r ist entweder ein array mit radii oder der equidistante
    % Abstand für die radii
    if length(delta_r) > 1
        r = delta_r;
    else
        Nr = rmax / delta_r;
        r = linspace(0, rmax, Nr+1)';  % Radial grid (column vector)
    end

    % t ist entweder ein array mit Zeiten, für die u als Matrix ausgegeben
    % werden soll (Format r x t), oder ein einzelner Zeitpunkt für den u 
    % als array ausgegeben werden soll
    if length(t) > 1
        u = zeros(length(r), length(t));
        for t_idx = 1:length(t)
            u_at_t = 1 * ((1 / pi) / (1 + 4*1*D*t(t_idx)))^(dim/2) .* exp(-(1 * r.^2) / (1 + 4*1*D*t(t_idx)) + gamma * t(t_idx));
            u(:, t_idx) = u_at_t;
        end
    else
        % solution
        u = 1 * ((1 / pi) / (1 + 4*1*D*t))^(dim/2) .* exp(-(1 * r.^2) / (1 + 4*1*D*t) + gamma * t);  % since r = sqrt(x^2 + y^2), I use r^2 = |x|^2
    end

end