% Solver for the 2D reaction-diffusion equation with exponential growth.
% Spatial grid and initial condition can be switched, read the comments.

function [u, r] = solution_RD_2D_exp(t_end, gamma, D)  %, rmax, delta_t, delta_r, epsilon)
    % variable parameters (can be used as input arguments as well)
    rmax = 10;
    delta_t = 0.00012817;
    delta_r = 0.015625;
    epsilon = 0.0781;

    % fixed parameters (dimension-dependent & should only be modified if 1D or 3D settings are computed)
    dim = 2;
    d = 4;  % for the Crank-Nicolson matrices
    int_factor = 2 * pi;  % for the integral

    % time steps
    nt = int32(t_end / delta_t);   % Number of time steps

    % spatial grid (column vector)
    % R_unif (base for the refinement, don't comment!)
    Nr_coarse = rmax / delta_r;
    r_coarse = linspace(0, rmax, Nr_coarse+1);
    r = r_coarse';  % R_unif (this line can be commented if a stretching is used)

    % % grid stretchings (if desired)
    % % R_1 and R_2 (comment the single line with the R_i you don't want)
    % r_coarse = r_coarse(r_coarse > 1);  % preparation
    % Nr_fine = ceil(1 / (delta_r / 10));  % preparation
    % r_fine = linspace(0, 1, Nr_fine + 1);  % R_1 
    % % r_fine = r_fine.^3;  % R_2
    % r = [r_fine, r_coarse]';  % finished grid
    % 
    % % R_3
    % r = (1/rmax * r_coarse.^2)';
    % 
    % % R_4
    % r = (1/(rmax^2) * r_coarse.^3)';
    
    % compute the potentially non-uniform distances between the grid points
    delta_r = diff(r);
    Nr = length(delta_r);

    nt = int32(t_end / delta_t);  % Number of time steps

    % % Initial condition
    % % Dirac-Delta approximation i.c. (same order as in the publication)
    % % un-comment the desired one
    % u = usual_u0(epsilon, r, dim, delta_r);
    % u = fraction_u0(epsilon, r, dim, delta_r);
    % u = rectangle_u0(epsilon, r, dim, delta_r);
    % u = power_u0(epsilon, r, dim, delta_r);
    % % Gaussian curve i.c.
    u = u0_gaussian(epsilon, r, dim, delta_r);

    % check if the initial condition has the correct mass (trapezoidal
    % rule)
    init_mass = int_factor * trapz(r, r.^(dim-1) .* u);
    integral_i_c = init_mass;

    % save all time steps
    u_all = zeros(Nr+1, nt);

    % Prepare the Crank-Nicolson method
    A = zeros(Nr+1, Nr+1);  % Matrix for the implicit method
    B = zeros(Nr+1, Nr+1);  % Matrix for the explicit method

    for i = 2:Nr
        dri = delta_r(i-1) + delta_r(i);  % will be used a lot, distance between r_{i-1} and r_{i+1}
        A(i,i-1) = -D * (1/(dri * delta_r(i-1)) - d/(2*r(i)*dri));        
        A(i,i)   = 1/delta_t + D/dri * (1/delta_r(i) + 1/delta_r(i-1)) - gamma/2;
        A(i,i+1) = -D * (1/(dri * delta_r(i)) + d/(2*r(i)*dri)); 

        B(i,i-1) = D * (1/(dri * delta_r(i-1)) - d/(2*r(i)*dri));
        B(i,i)   = 1/delta_t - D/dri * (1/delta_r(i) + 1/delta_r(i-1)) + gamma/2;
        B(i,i+1) = D * (1/(dri * delta_r(i)) + d/(2*r(i)*dri)); 
    end

    % Boundary conditions (Neumann BC at r=0, Dirichlet at r=rmax)
    A(1,1) = 1/delta_t + D / delta_r(1)^2; A(1,2) = - D / delta_r(1)^2;  % Symmetry condition at r=0
    B(1,1) = 1/delta_t - D / delta_r(1)^2; B(1,2) = D / delta_r(1)^2;  % Symmetry condition at r=0

    A(end,end) = 1; A(end,end-1) = 0;  % Dirichlet condition at r=rmax
    B(end,end) = 1; B(end,end-1) = 0;  % Dirichlet condition at r=rmax

    % Time-stepping loop
    for n = 1:nt
        rhs = B * u;
        u = A \ rhs;
        
        u_all(:, n) = u;
    end

    % % integral over the final solution to check mass development
    % integral_tend = int_factor * trapz(r, r.^(dim-1) .* u);
    % final_mass = integral_tend;
    % if length(r(r >= 1 & r <= 2)) > 1  % stellt sicher dass im gewählten Intervall mehr als ein grid point ist
    %     r_small_range = r(r >= 1 & r <= 2);
    %     int_u_low_r = int_factor * trapz(r_small_range, r_small_range.^(dim-1) .* u(r >= 1 & r <= 2));
    % else
    %     int_u_low_r = 0;
    % end
    % 
    % disp(['The numerical solution has an initial mass: ', num2str(integral_i_c), ', final mass: ', num2str(integral_tend), ', final mass for 1<r<2: ', num2str(int_u_low_r)])

    function u0 = usual_u0(epsilon, r, dim, delta_r)
        u0 = zeros(length(r), 1);
        
        % OPTION 1 
        % rectangular function
        % r_less_or_slightly_larger_than_eps = r(1:find(r > epsilon, 1));
        % % solves 1 = u0 * \sum (r_{i+1}^2 - r_i^2)
        % sum = r_less_or_slightly_larger_than_eps(end)^dim - r_less_or_slightly_larger_than_eps(1)^dim;  % many things cancel in the sum over the r_i
        % u0(1:find(r > epsilon, 1)) = 1 / sum;
    
        % OPTION 2
        % rectangular * n-dimensional delta 
        u0(r <= epsilon) = 1 / epsilon;
    
        % use the n-dimensional delta-function 
        if dim == 2
            u0(2:end) = 1./(2*pi*r(2:end)).*u0(2:end);
        elseif dim == 3
            u0(2:end) = 1./(4*pi*r(2:end).^2).*u0(2:end);
        end
    
        % integral over the "valid" part of the i.c. (i.e. u(2:end))
        if dim == 2
            int_factor = 2 * pi;  % for the integral
        elseif dim == 3
            int_factor = 4 * pi;
        end
    
        % trapezoidal rule
        int_ohne_r1_und_r2 = int_factor * trapz(r(3:end), r(3:end).^(dim-1) .* u0(3:end));
    
        rest = 1 - int_ohne_r1_und_r2;
        if rest < 0
            stack = dbstack;
            disp(['Problem in Function ', stack(1).name, ', epsilon = ', num2str(epsilon)]);
            disp('-> integral over the i.c. already > 1, even though u0(1), u0(2) have not been considered yet.')
            % factor = 1.1;
            % while rest < 0
            %     disp(['u0 wurde neu gestartet mit y-intercept 1 / (' num2str(factor) ' * epsilon).'])
            %     % nochmal mit niedrigerer Rechteck-delta approx versuchen
            %     u0(r <= epsilon) = 1 / (factor * epsilon);
            %     % use the n-dimensional delta-function 
            %     if dim == 2
            %         u0(2:end) = 1./(2*pi*r(2:end)).*u0(2:end);
            %     elseif dim == 3
            %         u0(2:end) = 1./(4*pi*r(2:end).^2).*u0(2:end);
            %     end
            % 
            %     % integral over the "valid" part of the i.c. (i.e. u(2:end))
            %     if dim == 2
            %         int_factor = 2 * pi;  % for the integral
            %     elseif dim == 3
            %         int_factor = 4 * pi;
            %     end
            % 
            %     % trapezoidal rule
            %     int_ohne_r1_und_r2 = int_factor * trapz(r(3:end), r(3:end).^(dim-1) .* u0(3:end));
            %     rest = 1 - int_ohne_r1_und_r2;
            %     factor = factor + 0.1;
            % end
        end
    
        % wenn es keine Probleme (mehr) mit rest gibt
        % trapezoidal rule
        % for uniform grids
        %alpha = (rest/(int_factor * 0.5 * delta_r) - u0(3) * r(3)^(dim-1)) / (r(1)^(dim-1) + 2 * r(2)^(dim-1));
        % for non-uniform grids
        alpha = (rest/(0.5 * int_factor) - delta_r(2) * u0(3) * r(3)^(dim-1)) / (delta_r(1) * (r(1)^(dim-1) + r(2)^(dim-1)) + delta_r(2) * r(2)^(dim-1)); 
        u0(1) = alpha;
        u0(2) = alpha;
        
        % OPTION 3
        % triangle function (height 2 / epsilon, width epsilon => 1D-area = 1)
        % r_less_than_eps = r(r <= epsilon);
        % for i = 1:length(r_less_than_eps)
        %     u0(i) = 2 / epsilon - 2 / (epsilon^2) * r(i);
        % end
    end
    
    
    function u0 = fraction_u0(epsilon, r, dim, delta_r)
        % 1/pi * lim_{eps -> 0} eps/(r^2 + eps^2)
        u0 = ones(length(r), 1);
    
        % use the n-dimensional delta-function 
        if dim == 2
            u0(2:end) = 1/pi * 1./(2*pi*r(2:end)).* epsilon./(r(2:end).^2 + epsilon^2) .* u0(2:end);
        elseif dim == 3
            u0(2:end) = 1/pi * 1./(4*pi*r(2:end).^2).* epsilon./(r(2:end).^2 + epsilon^2) .* u0(2:end);
        end
    
        % um das Risiko von "Masse aus dem Nichts" zu vermeiden, ab Rmax/2
        % einfach alles gleich 0 setzen
        half = ceil(length(r)/2);
        u0(half:end) = 0;
    
        % integral over the "valid" part of the i.c. (i.e. u(2:end))
        if dim == 2
            int_factor = 2 * pi;  % for the integral
        elseif dim == 3
            int_factor = 4 * pi;
        end
    
        % trapezoidal rule
        int_ohne_r1_und_r2 = int_factor * trapz(r(3:end), r(3:end).^(dim-1) .* u0(3:end));
    
        rest = 1 - int_ohne_r1_und_r2;
        if rest < 0
            stack = dbstack;
            disp(['Problem in Function ', stack(1).name, ', epsilon = ', num2str(epsilon)]);
            disp('-> integral over the i.c. already > 1, even though u0(1), u0(2) have not been considered yet.')
        end
    
        % wenn es keine Probleme (mehr) mit rest gibt
        % trapezoidal rule
        % for uniform grids
        %alpha = (rest/(int_factor * 0.5 * delta_r) - u0(3) * r(3)^(dim-1)) / (r(1)^(dim-1) + 2 * r(2)^(dim-1));
        % for non-uniform grids
        alpha = (rest/(0.5 * int_factor) - delta_r(2) * u0(3) * r(3)^(dim-1)) / (delta_r(1) * (r(1)^(dim-1) + r(2)^(dim-1)) + delta_r(2) * r(2)^(dim-1)); 
        u0(1) = alpha;
        u0(2) = alpha;
    end
    

    function u0 = rectangle_u0(epsilon, r, dim, ~)
        % eigentlich ein Zylinder mit Radius 1/2 epsilon und einer Höhe die so
        % passt dass das Volumenintegral gleich 1 ist
        u0 = zeros(length(r), 1);
    
        % integrating factor for later
        if dim == 2
            int_factor = 2 * pi;  % for the integral
        elseif dim == 3
            int_factor = 4 * pi;
        end
    
        % last radius in r that is <= epsilon
        idx = find(r <= 0.5*epsilon, 1, 'last');
    
        % bei sehr grobem Gitter sicherstellen dass wenigstens zwei
        % Gitterpunkte gleich 1/epsilon_star sein können
        if isempty(idx) || idx == 1
            idx = 2;
        end
        
        % distances in r
        diff_r = diff(r);
    
        % height of the "rectangle" based on the trapez rule with volume
        % element and a non-uniform grid
        epsilon_star = 0.5 * int_factor * (trapz(r(1:idx), r(1:idx).^(dim-1)) + diff_r(idx) * r(idx)^(dim-1));
    
        u0(1:idx) = 1 / epsilon_star;
    end


    function u0 = power_u0(epsilon, r, dim, delta_r)
        % lim_{eps -> 0} eps * r^(eps-1)
        u0 = ones(length(r), 1);
    
        % use the n-dimensional delta-function 
        if dim == 2
            u0(2:end) = epsilon * 1./(2*pi*r(2:end)).* r(2:end).^(epsilon - 1) .* u0(2:end);
        elseif dim == 3
            u0(2:end) = epsilon * 1./(4*pi*r(2:end).^2).* r(2:end).^(epsilon - 1) .* u0(2:end);
        end
    
        % um das Risiko von "Masse aus dem Nichts" zu vermeiden, ab Rmax/2
        % einfach alles gleich 0 setzen
        half = ceil(length(r)/2);
        u0(half:end) = 0;
    
        % integral over the "valid" part of the i.c. (i.e. u(2:end))
        if dim == 2
            int_factor = 2 * pi;  % for the integral
        elseif dim == 3
            int_factor = 4 * pi;
        end
    
        % trapezoidal rule
        int_ohne_r1_und_r2 = int_factor * trapz(r(3:end), r(3:end).^(dim-1) .* u0(3:end));
    
        rest = 1 - int_ohne_r1_und_r2;
        if rest < 0
            stack = dbstack;
            disp(['Problem in Function ', stack(1).name, ', epsilon = ', num2str(epsilon)]);
            disp('-> integral over the i.c. already > 1, even though u0(1), u0(2) have not been considered yet.')
        end
    
        % wenn es keine Probleme (mehr) mit rest gibt
        % trapezoidal rule
        % for uniform grids
        %alpha = (rest/(int_factor * 0.5 * delta_r) - u0(3) * r(3)^(dim-1)) / (r(1)^(dim-1) + 2 * r(2)^(dim-1));
        % for non-uniform grids
        alpha = (rest/(0.5 * int_factor) - delta_r(2) * u0(3) * r(3)^(dim-1)) / (delta_r(1) * (r(1)^(dim-1) + r(2)^(dim-1)) + delta_r(2) * r(2)^(dim-1)); 
        u0(1) = alpha;
        u0(2) = alpha;
    end
    
    
    function u0 = u0_gaussian(~, r, dim, ~)
        % I.c. passend zur analytical solution, a = 1, rho0 = 1.
        u0 = (1/pi)^(dim/2) * exp(-r.^2);
    end
end