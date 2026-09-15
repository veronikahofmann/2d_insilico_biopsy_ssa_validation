% Script to simulate the biopsies based on solution_RD_2D_log.m and the 
% analytical solution of the 2D RD-equation with exponential growth and 
% Dirac-Delta initial condition. Before running the script make sure to 
% have grid "R_1" and initial condition "usual_u0" set in the solver for 
% best performance.
% If you want to save the cell positions as excel files, indicate this
% below.

save_coordinates = true;
version = 00;
folder = 'insilico_biopsies_from_that_weird_paper';

% parameter combinations (order: gamma, D)
V = [0.3, 5000];  % minimum proliferation, minimum diffusion
W = [0.7, 150000];  % max prolif, max diff
X = [0.7, 5000];  % max prolif, min diff
Y = [0.3, 150000];  % min prolif, max diff

% solver parameters
maxT = 6;
maxR = 10;
dim = 2;
delta_r = 0.015625;

% random number generator
seed = 0;
rng(seed,"twister") 

% physical data & parameters
% cell radius can be given here as a scalar or later, once r is known, 
% as a spatial distribution from 0 to maxR
cell_radius = 11.0;  % in micrometers
length_unit = 1000;  % in micrometers, unit in which the calculations are performed

% noise parameter (of the cell size) -> has to be log-scale if log-normal
% distributed noise is selected further below!!
sd_size = 0.2;  % same unit as cell size, i.e. micrometers

% noise parameters (of the cell position) 
sd_noise = 0.02;  % noise standard deviation, unit: length units
sd_noise_microm = sd_noise*length_unit;  % sd_noise in micrometers
sd_count = ceil(3*sd_noise_microm/(delta_r*length_unit));  % count how many d_r's are covered by 3 stds

% data management
letters = ["V", "W", "X", "Y"];
parameter_collection = dictionary(letters, {V, W, X, Y});

if ~exist(folder, 'dir')
   mkdir(folder)
end

for ind_letter = 1:length(letters)
    % load parameter combination together with parameter combination
    % identifier (= letter)
    letter = letters(ind_letter);  
    parameters_cellobj = parameter_collection(letter);
    parameters = parameters_cellobj{1};
    gamma = parameters(1);
    D = parameters(2)/(length_unit^2);  % conversion to mm²/day

    disp(strcat('Biopsy for Parameter Combination ', " ", letter, ' is about to be computed...'));

    % file management in case saving is desired
    letterChar = char(letter);
    filename_xy = strcat(folder, '/group_', letterChar, '.xlsx');

    % compute the solution of the RD-equation -> choose the desired option
    % normalized exponential growth
    [u_unnormed, r] = analytical_solution_RD_exp_Dirac(dim, maxT, maxR, gamma, D, delta_r);  
    [u_log_for_norm, ~] = solution_RD_2D_log(maxT, gamma, D);  % compute the solution of the logistic-growth RD-equation for the "norm"
    u_norm = max(u_unnormed) / max(u_log_for_norm);
    u = u_unnormed / u_norm;
    % % logistic growth
    % [u_log_for_norm, r] = solution_RD_2D_log(maxT, gamma, D);

    % % CELL RADIUS VARIABLE UNDER r
    % % now that we have the array r, we can define the cell radius as a
    % % distribution along r (e.g. if we expect cells close to the center of
    % % the tumor to have a different mean size than further away)
    % % (the scaling works different for each type of distribution, for some
    % % distributions you need a max. value, for others a mean)
    % lambda = 0.1; 
    % scaling_factor = 10 / lambda;
    % cell_radius = scaling_factor * flip(exppdf(r, 1/lambda));

    % ALTERNATIVE: CELL RADIUS INDEPENDENT OF r
    cell_radius = cell_radius;

    % plot the cell size distribution to make sure that it is set up
    % correctly
    
    % % noise-band if Gaussian (width = 2 stds)
    % if isscalar(cell_radius)
    %     upper = ones(size(r)) * cell_radius + 2 * sd_size;
    %     lower = ones(size(r)) * cell_radius - 2 * sd_size;
    % else
    %     upper = cell_radius + 2 * sd_size;
    %     lower = cell_radius - 2 * sd_size;
    % end

    % noise band if log-normal (width = 95% bc. 2 std does not make sense due to assymmetry)
    crit_val = norminv(1 - 0.05/2);  % skaliert das Konfidenzintervall: wenn man ein 1 - alpha Konfidenzintervall möchte, muss man hier 1 - alpha/2 schreiben
    % log-normal dist., also ln(cell_radius) \sim N(mu, sigma²), hat andere
    % parameter in der log-domain (mu, sigma) als in der linearen "Daten-domain" 
    sigma_sq = sd_size^2;  %log(1 + (sd_size ./ cell_radius).^2); ACHTUNG STD IST BEI VIELEN PAPERS SCHON LOG SCALE
    mu = log(cell_radius) - sigma_sq ./ 2;
    if isscalar(cell_radius)
        upper = exp(ones(size(r)) * mu + sqrt(sigma_sq) * crit_val);
        lower = exp(ones(size(r)) * mu - sqrt(sigma_sq) * crit_val);
    else
        upper = exp(mu + sqrt(sigma_sq) * crit_val);
        lower = exp(mu - sqrt(sigma_sq) * crit_val);
    end

    a = [r', fliplr(r')];
    b = [upper', fliplr(lower')];

    figure;
    fill(a, b, [0.9 0.9 0.9], 'DisplayName', 'noise (95% confidence)', 'EdgeColor', 'none', 'FaceAlpha', 0.5);  % 'noise (2 standard deviations)'
    hold on
    if isscalar(cell_radius)
        plot(r, ones(size(r)) * cell_radius, 'LineWidth', 2, 'DisplayName', 'mean cell radius')
    else
        plot(r, cell_radius, 'LineWidth', 2, 'DisplayName', 'mean cell radius')
    end
    xlabel('domain radius $r$ in $1000 \mu$m', 'Interpreter', 'latex');
    ylabel('cell radius in $\mu$m', 'Interpreter', 'latex');
    title('Cell Size Distribution in the Spatial Domain');
    legend();
    ylim([0 15])
    hold off
    drawnow;

    % compute theoretical maximum number of cells in the biopsy
    num_points = AvailableSpace(r, length_unit, cell_radius);

    % Create arrays to store 2D coordinates and concentrations
    x_slice = [];
    y_slice = [];
    sizes_slice = [];
    u_values = [];
    alpha_values = [];

    % Loop through each radial distance and generate random points
    Nr = length(r);

    for i = 1:Nr
        % disp(['Computing point cohort ', num2str(i), '/', num2str(Nr)])

        % Generate random angle for the polar coordinates
        phi = 2*pi*rand(num_points(i), 1);
        % Generate random radius
        radius = r(i) .* (1 + sd_noise .* randn(num_points(i), 1));
        % Convert spherical to Cartesian coordinates
        x = radius .* cos(phi);
        y = radius .* sin(phi);

        % optional: add log-normal noise to the cell radius
        % distribution in cell_radius + log-normally distr. noise
        if isscalar(cell_radius)
            sizes = lognrnd(mu, sqrt(sigma_sq), size(y));
        else
            sizes = lognrnd(mu(i), sqrt(sigma_sq), size(y)); 
        end

        % start collision avoidance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % % test if the point lies to close to another point (cells cannot
        % % overlap)
        % % distances between the points in this cohort
        % dx = x' - x;
        % dy = y' - y;
        % distance = sqrt(dx.^2 + dy.^2);  % pairwise distance matrix
        % % % früher (immer gleiche Zellgröße)
        % % point_invalid = distance < 2*cell_radius/length_unit;  % only continue if no intersection
        % % jetzt (immer mindestens mit bissl noise, sodass alle Zellen unterschiedlich groß sind)
        % combined_cell_rad = (sizes + sizes')./length_unit;
        % point_invalid = distance < combined_cell_rad;
        % point_invalid = point_invalid - diag(diag(point_invalid));  % Diagonale wird weggelassen, da eine Zelle von sich selbst natürlich 0 Entfernung hat
        % points_to_delete = false(1, length(x));  % Initialize delete flag array
        % 
        % for j = 1:length(x)
        %     % Find indices j > i where distances(i,j) is below threshold
        %     close_points = find(point_invalid(j, j+1:end)) + j;
        %     % Mark the first point in each close pair (only once)
        %     if ~isempty(close_points)
        %         points_to_delete(j) = true;  % Mark point i for deletion
        %         point_invalid(j, :) = 0;       % Clear row i to avoid redundant deletions
        %         point_invalid(:, j) = 0;       % Clear column i to avoid cascading deletions
        %     end
        % end
        % x = x(~points_to_delete);
        % y = y(~points_to_delete);
        % sizes = sizes(~points_to_delete);
        % 
        % % distance to previously created points
        % if ~isempty(x_slice)  % if the list of previously created pts exists
        %     % point_valid = false(num_points, 1);
        %     % die anzahl der zu testenden Pkte hängt von der std der
        %     % Positionierung und der maximalen Zellgröße ab
        %     max_cell_radius = max(cell_radius) + 3 * sd_size;  % [µm]
        %     % die max. # der Kreise die vom Radius einer maximal großen
        %     % Zelle besetzt sein könnten
        %     max_nr_kreise_cellsize = ceil(max_cell_radius / (delta_r*length_unit));  % [µm]
        %     nr_pts_to_be_tested = min(length(x_slice), 2*(sd_count + max_nr_kreise_cellsize)*num_points(i));  % siehe sticky note April 2025
        %     min_idx_tbt = length(x_slice) - nr_pts_to_be_tested + 1;
        %     indices_to_be_tested = min_idx_tbt:length(x_slice);
        %     for j = indices_to_be_tested
        %         distance = sqrt((x_slice(j) - x).^2 + (y_slice(j) - y).^2);
        %         % % früher (immer gleiche Zellgröße)
        %         % point_valid = distance >= 2*cell_radius/length_unit;  % only continue if no intersection
        %         % jetzt (immer mindestens mit bissl noise, sodass alle Zellen unterschiedlich groß sind)
        %         combined_cell_rad = (sizes_slice(j) * ones(size(sizes)) + sizes)./length_unit;
        %         point_valid = distance >= combined_cell_rad;
        %         x = x(point_valid); % get rid of invalid points immediately
        %         y = y(point_valid);
        %         sizes = sizes(point_valid);
        %     end
        % else  % if it's the first time we create pts
        %     point_valid = true(length(x), 1);
        % end

        % end collision avoidance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Store points and corresponding concentrations
        for pt = 1:length(x)
            rnd = rand();  % random number between 0 and 1, uniform distribution
            if rnd < u(i)  % only create a point by chance
                x_slice = [x_slice; x(pt)];
                y_slice = [y_slice; y(pt)];
                sizes_slice = [sizes_slice; sizes(pt)];
                u_values = [u_values; u(i)];

                % Set transparency level based on u(r): 1 when u=1, and 0 when u=0
                alpha_values = [alpha_values; u(i)];
            end
        end
    end

    % retrieving the concentration from the point pattern
    dists_to_zero = sqrt(x_slice.^2 + y_slice.^2);

    fractions = zeros(size(r));

    for ind = 1:length(r)
        ri = r(ind);
        if ri == 0
            area = pi * (delta_r/2)^2;  % area of the partial domain
            points_in_range = (dists_to_zero < delta_r/2);  % cells with nucleus inside partial domain
        else
            area = pi * ((ri + delta_r/2)^2 - (ri - delta_r/2)^2);
            points_in_range = (dists_to_zero >= ri - delta_r/2) & (dists_to_zero < ri + delta_r/2);
        end
        % get the positions and sizes of the cells that are in the current ring
        x_in_range = x_slice(points_in_range);
        y_in_range = y_slice(points_in_range);
        sizes_in_range = sizes_slice(points_in_range) ./ length_unit;

        if isempty(x_in_range)
            area_cells = 0;
        else
            area_cells = pi * sum(sizes_in_range.^2); % Gesamtfläche der Zellen berechnen, deren Nuclei im aktuellen Ring sind (egal ob was aus dem Ring rauslappt oder nicht)
        end
        
        fractions(ind) = area_cells / area;  % sizes_in_range µm, r in LE
    end

    % save for a plot with the results from multiple parameter combinations
    areafractions{ind_letter} = fractions;
    rs{ind_letter} = r;
    us{ind_letter} = u;
    gammas(ind_letter) = gamma;
    Ds(ind_letter) = D;
    maxts(ind_letter) = maxT;

    % plot the resulting concentration and biopsy
    figure;
    tiles = tiledlayout(1, 2);
    nexttile
    plot(r, u_unnormed, 'LineWidth', 2, 'LineStyle', '--')
    hold on
    plot(r, u, 'LineWidth', 2, 'LineStyle', ':');
    plot(r, fractions, 'LineWidth', 1)
    legend({'$u_h$', '$\| u_h \|$', 'in-silico biopsy'}, 'Interpreter', 'latex')
    title('Solution of the 2D RD-equation', 'Interpreter', 'latex')
    xlabel("$r$ ($\times 1000 \: \mu$m)", 'Interpreter', 'latex'); 
    ylabel("cell density $u(r, t)$", 'Interpreter', 'latex');
    ylim([0 1])
    grid on
    hold off
    nexttile
    u_values(u_values > 1) = 1;  % im Code wird u > 1 ja sowieso wie u = 1 behandelt, nur für die colorbar ist das nicht selbstverständlich
    scatter(x_slice, y_slice, 10, u_values, 'filled');
    grid on
    colorbar
    clim([0, 1]);
    title('In-silico biopsy generated from $\| u_h \|$', 'Interpreter', 'latex')

    %titleString = sprintf('Parameter Combination %s: $t = %.1f$, $\gamma = %.2f$, $D = %.4f$', letter, maxT, gamma, D);
    titleString = ['Parameter Combination v', num2str(version), '/', letterChar, ': $\gamma = ', num2str(gamma), '$, $D = ', num2str(D), '$ ($t = ', num2str(maxT), '$)'];
    title(tiles, titleString, 'Interpreter', 'latex');

    drawnow
    
    % save the slice
    if save_coordinates
        t_array = repmat('t', length(x_slice), 1);
        T = table(t_array, x_slice*length_unit, y_slice*length_unit, 'VariableNames', ["class", "x", "y"]);
        writetable(T, filename_xy, 'FileType','spreadsheet');
    end  

end

% Plot the deviations of the measured concentrations from the (normalized) u
figure;
titles = [];
tiles = tiledlayout(1, length(letterlist));
for i = 1:length(letterlist)
    nexttile;
    plot(rs{i}, us{i}, 'LineWidth', 2);
    hold on
    plot(rs{i}, areafractions{i}, 'LineWidth', 2, 'LineStyle', ':')
    legend({'computed density $\| u_h \|$', 'measured density in-silico biopsy'}, 'Interpreter', 'latex')
    title('Solution of the 2D RD-equation', 'Interpreter', 'latex')
    xlabel("$r$ ($\times 1000 \: \mu$m)", 'Interpreter', 'latex'); 
    ylabel("cell density $u(r, t)$", 'Interpreter', 'latex');
    titlestr = ['$\gamma = ', num2str(gammas(i)), '$, $D = ', num2str(Ds(i)), '$ ($t = ', num2str(maxts(i)), '$)'];
    title(titlestr, 'Interpreter', 'latex')
    ylim([0 1])
    xlim([0 6])
    grid on
    hold off
end

% helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N = AvailableSpace(r, length_unit, cell_radius)
% calculates how many cells fit on the circle with the current radius
    r_in_mu = length_unit*r;  % r in the physical unit (µm)
    circumference = 2 * r_in_mu * pi;
    if isscalar(cell_radius)  % if cell size is a const. value 
        N = ceil(circumference/(2*cell_radius));
    else  % if cell size is a size distribution over r
        N = ceil(circumference./(2*cell_radius));
    end
    N(1) = 1;
end

