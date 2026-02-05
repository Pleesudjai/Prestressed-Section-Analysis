function plotMultipleStressStrain(results, x_locations, options)
% PLOTMULTIPLESTRESSSTRAIN - Plot stress/strain distributions at multiple locations
%
% Compares stress and strain distributions at several locations along the beam
%
% Syntax:
%   plotMultipleStressStrain(results, x_locations)
%   plotMultipleStressStrain(results, x_locations, options)
%
% Inputs:
%   results     - Analysis resultsbbb structure from analyzePrestressedBeam
%   x_locations - Array of locations (in) or fractions of span (0-1)
%                 Default: [0.1, 0.25, 0.5, 0.75, 0.9]
%   options     - Optional structure with fields:
%                 .plot_type - 'stress', 'strain', or 'both' (default: 'stress')
%                 .num_points - Number of points through depth (default: 50)
%
% Example:
%   results = analyzePrestressedBeam(beam, section, materials, prestress, reinforcement, loads);
%   plotMultipleStressStrain(results, [0.1, 0.25, 0.5]);
%   plotMultipleStressStrain(results, [60, 150, 300, 450, 540]);

%% Parse inputs
if nargin < 2 || isempty(x_locations)
    x_locations = [0.1, 0.25, 0.5, 0.75, 0.9];
end

if nargin < 3
    options = struct();
end

if ~isfield(options, 'plot_type'), options.plot_type = 'stress'; end
if ~isfield(options, 'num_points'), options.num_points = 50; end

%% Setup
L = results.L;
x = results.x;
section = results.section;
materials = results.materials;

yb = section.yb;
yt = section.yt;
A = section.A;
Ix = section.Ix;
Ec = materials.Ec;

% Convert x_locations to actual positions
x_actual = zeros(size(x_locations));
for i = 1:length(x_locations)
    if x_locations(i) >= 0 && x_locations(i) <= 1
        x_actual(i) = x_locations(i) * L;
    else
        x_actual(i) = x_locations(i);
    end
end

n_locations = length(x_actual);

%% Create y-array through depth
y_from_bottom = linspace(0, yb + yt, options.num_points);
y_from_centroid = y_from_bottom - yb;

%% Calculate stress/strain at each location
stress_all = zeros(options.num_points, n_locations);
strain_all = zeros(options.num_points, n_locations);
P_all = zeros(1, n_locations);
e_all = zeros(1, n_locations);
M_all = zeros(1, n_locations);

for i = 1:n_locations
    [~, idx] = min(abs(x - x_actual(i)));
    
    P = results.P(idx);
    e = results.e(idx);
    M = results.M(idx);
    
    P_all(i) = P;
    e_all(i) = e;
    M_all(i) = M;
    
    % Total stress at each depth
    f_prestress = -P/A + P*e*y_from_centroid/Ix;
    f_external = -M*y_from_centroid/Ix;
    f_total = f_prestress + f_external;
    
    stress_all(:, i) = f_total;
    strain_all(:, i) = f_total / Ec * 1e6;  % microstrain
end

%% Create figure
figure('Name', 'Multi-Location Stress-Strain', 'Position', [100, 100, 1400, 600], 'Color', 'w');

colors = lines(n_locations);

%% Plot stress
if strcmp(options.plot_type, 'stress') || strcmp(options.plot_type, 'both')
    if strcmp(options.plot_type, 'both')
        subplot(1, 2, 1);
    end
    
    hold on;
    grid on;
    box on;
    
    h_section = max(y_from_bottom);
    
    % Zero line
    plot([0, 0], [0, h_section], 'k-', 'LineWidth', 1);
    
    % Centroid line
    plot(xlim, [yb, yb], 'k--', 'LineWidth', 0.5);
    
    % Plot each location
    legend_entries = cell(1, n_locations);
    h_lines = zeros(1, n_locations);
    for i = 1:n_locations
        h_lines(i) = plot(stress_all(:, i), y_from_bottom, '-', ...
            'Color', colors(i,:), 'LineWidth', 2);
        legend_entries{i} = sprintf('x = %.0f in (%.0f%%)', x_actual(i), x_actual(i)/L*100);
    end
    
    % Allowable limits
    stresses = results.stresses;
    fc_comp = stresses.fc_allow_compression;
    fc_tens = stresses.fc_allow_tension;
    
    plot([fc_comp, fc_comp], [0, h_section], 'k:', 'LineWidth', 1.5);
    plot([fc_tens, fc_tens], [0, h_section], 'k:', 'LineWidth', 1.5);
    
    xlabel('Stress, f (ksi)', 'FontSize', 11);
    ylabel('Distance from bottom (in)', 'FontSize', 11);
    title('Stress Distribution at Multiple Locations', 'FontWeight', 'bold', 'FontSize', 12);
    legend(h_lines, legend_entries, 'Location', 'best', 'FontSize', 9);
    
    % Adjust limits
    curr_xlim = xlim;
    xlim([min(curr_xlim(1)*1.1, fc_comp*1.1), curr_xlim(2)*1.1 + 0.2]);
end

%% Plot strain
if strcmp(options.plot_type, 'strain') || strcmp(options.plot_type, 'both')
    if strcmp(options.plot_type, 'both')
        subplot(1, 2, 2);
    end
    
    hold on;
    grid on;
    box on;
    
    h_section = max(y_from_bottom);
    
    % Zero line
    plot([0, 0], [0, h_section], 'k-', 'LineWidth', 1);
    
    % Centroid line
    plot(xlim, [yb, yb], 'k--', 'LineWidth', 0.5);
    
    % Plot each location
    legend_entries = cell(1, n_locations);
    h_lines = zeros(1, n_locations);
    for i = 1:n_locations
        h_lines(i) = plot(strain_all(:, i), y_from_bottom, '-', ...
            'Color', colors(i,:), 'LineWidth', 2);
        legend_entries{i} = sprintf('x = %.0f in (%.0f%%)', x_actual(i), x_actual(i)/L*100);
    end
    
    % Crushing strain limit
    eps_crush = -3000;
    plot([eps_crush, eps_crush], [0, h_section], 'k:', 'LineWidth', 1.5);
    
    xlabel('Strain, ε (microstrain)', 'FontSize', 11);
    ylabel('Distance from bottom (in)', 'FontSize', 11);
    title('Strain Distribution at Multiple Locations', 'FontWeight', 'bold', 'FontSize', 12);
    legend(h_lines, legend_entries, 'Location', 'best', 'FontSize', 9);
end

%% Super title
sgtitle('Stress and Strain Distribution Comparison', 'FontSize', 14, 'FontWeight', 'bold');

%% Print summary table
fprintf('\n============================================================\n');
fprintf('  STRESS-STRAIN SUMMARY AT MULTIPLE LOCATIONS\n');
fprintf('============================================================\n');
fprintf('%-10s %-10s %-10s %-10s %-12s %-12s\n', ...
    'x (in)', 'x/L (%)', 'P (kips)', 'e (in)', 'f_top (ksi)', 'f_bot (ksi)');
fprintf('------------------------------------------------------------\n');
for i = 1:n_locations
    fprintf('%-10.1f %-10.1f %-10.1f %-10.2f %-12.4f %-12.4f\n', ...
        x_actual(i), x_actual(i)/L*100, P_all(i), e_all(i), ...
        stress_all(end, i), stress_all(1, i));
end
fprintf('============================================================\n\n');

end
