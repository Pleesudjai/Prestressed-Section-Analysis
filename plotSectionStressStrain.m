function plotSectionStressStrain(results, x_location, options)
% PLOTSECTIONSTRESSSTRAIN - Plot stress and strain distribution along section height
%
% Plots the variation of stress and strain through the cross-section depth
% at a specified location along the beam length.
%
% Prestress effect is separated into:
%   - Axial component:   f_axial   = -P/A  (uniform compression)
%   - Bending component: f_bending = P*e*y/I (linear, due to eccentricity)
%
% Syntax:
%   plotSectionStressStrain(results, x_location)
%   plotSectionStressStrain(results, x_location, options)
%
% Inputs:
%   results    - Analysis results structure from analyzePrestressedBeam
%   x_location - Location along beam (in) or fraction of span (0-1)
%   options    - Optional structure with fields:
%                .plot_type    - 'both', 'stress', 'strain' (default: 'both')
%                .show_limits  - Show allowable stress limits (default: true)
%                .num_points   - Number of points through depth (default: 50)
%                .show_components - Show prestress/external separately (default: true)
%                .figure_handle - Existing figure handle (optional)
%
% Example:
%   results = analyzePrestressedBeam(beam, section, materials, prestress, reinforcement, loads);
%   plotSectionStressStrain(results, 0.5);  % At midspan
%   plotSectionStressStrain(results, 300);  % At x = 300 in

%% Parse inputs
if nargin < 3
    options = struct();
end

% Default options
if ~isfield(options, 'plot_type'), options.plot_type = 'both'; end
if ~isfield(options, 'show_limits'), options.show_limits = true; end
if ~isfield(options, 'num_points'), options.num_points = 50; end
if ~isfield(options, 'show_components'), options.show_components = true; end

%% Determine x-location
L = results.L;
x = results.x;

% If x_location is between 0 and 1, treat as fraction of span
if x_location >= 0 && x_location <= 1
    x_loc = x_location * L;
else
    x_loc = x_location;
end

% Find nearest index in results arrays
[~, idx] = min(abs(x - x_loc));
x_actual = x(idx);

%% Extract section properties
section = results.section;
yb = section.yb;  % Distance from centroid to bottom
yt = section.yt;  % Distance from centroid to top
yc = section.yc;  % Centroid from bottom
A = section.A;
Ix = section.Ix;

%% Extract forces at this location
P = results.P(idx);           % Effective prestress force (kips)
e = results.e(idx);           % Eccentricity (in)
M = results.M(idx);           % External moment (kip-in)
N = results.N(idx);           % Axial force (kips)

%% Get material properties
materials = results.materials;
Ec = materials.Ec;            % Concrete modulus (ksi)

%% Calculate stress components at top and bottom fibers
% --- Prestress AXIAL component: f = -P/A (uniform compression) ---
f_axial = -P/A;  % Same at every fiber

f_axial_top = f_axial;
f_axial_bot = f_axial;

% --- Prestress BENDING component: f = P*e*y/I (linear) ---
f_prestress_bend_top = P*e*yt/Ix;
f_prestress_bend_bot = -P*e*yb/Ix;

% --- Combined prestress (for reference) ---
f_prestress_top = f_axial_top + f_prestress_bend_top;
f_prestress_bot = f_axial_bot + f_prestress_bend_bot;

% --- External load stresses: f = -M*y/I ---
f_external_top = -M*yt/Ix;
f_external_bot = M*yb/Ix;

% --- Total stresses ---
f_total_top = f_prestress_top + f_external_top;
f_total_bot = f_prestress_bot + f_external_bot;

%% Create arrays for plotting through depth
y_from_bottom = linspace(0, yb + yt, options.num_points);  % y from bottom of section
y_from_centroid = y_from_bottom - yb;  % y from centroid

% Stress distributions through depth
f_axial_arr       = -P/A * ones(size(y_from_centroid));       % Uniform axial
f_prebend_arr     = P*e*y_from_centroid/Ix;                   % Prestress bending
f_prestress_arr   = f_axial_arr + f_prebend_arr;              % Combined prestress
f_external_arr    = -M*y_from_centroid/Ix;                    % External load
f_total_arr       = f_prestress_arr + f_external_arr;         % Total

% Strain distributions (elastic)
epsilon_axial     = f_axial_arr / Ec;
epsilon_prebend   = f_prebend_arr / Ec;
epsilon_prestress = f_prestress_arr / Ec;
epsilon_external  = f_external_arr / Ec;
epsilon_total     = f_total_arr / Ec;

%% Find neutral axis location (where total stress = 0)
M_net = P*e - M;  % Net moment about centroid

if abs(M_net) > 1e-6
    y_na_from_centroid = (P/A) * Ix / M_net;
    y_na_from_bottom = yb + y_na_from_centroid;
    
    if y_na_from_bottom < 0 || y_na_from_bottom > (yb + yt)
        y_na_within_section = false;
    else
        y_na_within_section = true;
    end
else
    y_na_from_bottom = NaN;
    y_na_from_centroid = NaN;
    y_na_within_section = false;
end

% Calculate curvature
h_total = yb + yt;
curvature = (f_total_top - f_total_bot) / (Ec * h_total);  % 1/in

%% Create figure
if isfield(options, 'figure_handle') && ishandle(options.figure_handle)
    figure(options.figure_handle);
    clf;
else
    figure('Name', 'Section Stress-Strain Distribution', ...
        'Position', [100, 100, 1200, 600], 'Color', 'w');
end

%% Determine subplot layout
switch options.plot_type
    case 'both'
        num_cols = 2;
    otherwise
        num_cols = 1;
end

%% Plot stress distribution
if strcmp(options.plot_type, 'both') || strcmp(options.plot_type, 'stress')
    if strcmp(options.plot_type, 'both')
        subplot(1, 2, 1);
    end
    
    hold on;
    grid on;
    box on;
    
    h_section = max(y_from_bottom);
    
    % Plot stress components if requested
    if options.show_components
        % Prestress axial (uniform) - dashed green
        h1 = plot(f_axial_arr, y_from_bottom, 'g--', 'LineWidth', 1.5, ...
            'DisplayName', 'Prestress Axial (-P/A)');
        
        % Prestress bending (linear) - dash-dot green
        h2 = plot(f_prebend_arr, y_from_bottom, 'g-.', 'LineWidth', 1.5, ...
            'DisplayName', 'Prestress Bending (Pe \cdot y/I)');
        
        % Combined prestress - solid green
        h3 = plot(f_prestress_arr, y_from_bottom, 'g-', 'LineWidth', 2, ...
            'DisplayName', 'Prestress Combined');
        
        % External load - blue
        h4 = plot(f_external_arr, y_from_bottom, 'b-', 'LineWidth', 1.5, ...
            'DisplayName', 'External Load (-My/I)');
    end
    
    % Plot total stress - red, thickest
    h5 = plot(f_total_arr, y_from_bottom, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Total');
    
    % Fill compression and tension regions for total
    f_pos = max(f_total_arr, 0);
    f_neg = min(f_total_arr, 0);
    
    fill([f_pos, zeros(size(f_pos))], [y_from_bottom, fliplr(y_from_bottom)], ...
        [1, 0.7, 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([f_neg, zeros(size(f_neg))], [y_from_bottom, fliplr(y_from_bottom)], ...
        [0.7, 0.7, 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Zero line
    plot([0, 0], [0, h_section], 'k-', 'LineWidth', 1);
    
    % Centroid line
    plot(xlim, [yb, yb], 'k--', 'LineWidth', 0.8);
    text(max(xlim)*0.95, yb, ' C.G.', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'right', 'FontSize', 9);
    
    % Neutral axis
    if y_na_within_section
        plot(xlim, [y_na_from_bottom, y_na_from_bottom], 'm-', 'LineWidth', 2);
        text(max(xlim)*0.95, y_na_from_bottom, sprintf(' N.A. @ %.2f in', y_na_from_bottom), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
            'FontSize', 9, 'Color', 'm', 'FontWeight', 'bold');
    elseif ~isnan(y_na_from_bottom)
        if y_na_from_bottom < 0
            text(0.5, 0.05, sprintf('N.A. @ %.1f in (below section)', y_na_from_bottom), ...
                'Units', 'normalized', 'FontSize', 8, 'Color', 'm', ...
                'HorizontalAlignment', 'center');
        else
            text(0.5, 0.95, sprintf('N.A. @ %.1f in (above section)', y_na_from_bottom), ...
                'Units', 'normalized', 'FontSize', 8, 'Color', 'm', ...
                'HorizontalAlignment', 'center');
        end
    end
    
    % Allowable stress limits
    if options.show_limits
        stresses = results.stresses;
        fc_comp = stresses.fc_allow_compression;
        fc_tens = stresses.fc_allow_tension;
        
        plot([fc_comp, fc_comp], [0, h_section], 'k:', 'LineWidth', 1.5);
        text(fc_comp, h_section*0.02, sprintf(' f''_c = %.2f ksi', fc_comp), ...
            'VerticalAlignment', 'bottom', 'FontSize', 8, 'Rotation', 90);
        
        plot([fc_tens, fc_tens], [0, h_section], 'k:', 'LineWidth', 1.5);
        text(fc_tens, h_section*0.02, sprintf(' f_r = %.3f ksi', fc_tens), ...
            'VerticalAlignment', 'bottom', 'FontSize', 8, 'Rotation', 90);
    end
    
    % Annotations for extreme fiber stresses
    text(f_total_top, h_section, sprintf('  %.3f ksi', f_total_top), ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');
    text(f_total_bot, 0, sprintf('  %.3f ksi', f_total_bot), ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');
    
    % Labels
    xlabel('Stress, f (ksi)', 'FontSize', 11);
    ylabel('Distance from bottom (in)', 'FontSize', 11);
    title(sprintf('Stress Distribution at x = %.1f in (%.1f%% of span)', x_actual, x_actual/L*100), ...
        'FontWeight', 'bold', 'FontSize', 12);
    
    % Legend
    if options.show_components
        legend([h1, h2, h3, h4, h5], ...
            {'Prestress Axial (-P/A)', 'Prestress Bending (Pe \cdot y/I)', ...
             'Prestress Combined', 'External Load (-My/I)', 'Total'}, ...
            'Location', 'best', 'FontSize', 8);
    else
        legend(h5, {'Total'}, 'Location', 'best', 'FontSize', 9);
    end
    
    % Add text box with key values (now includes axial/bending breakdown)
    info_str = sprintf(['P = %.1f kips\n' ...
                       'e = %.2f in\n' ...
                       '-P/A = %.4f ksi\n' ...
                       'Pe = %.0f kip-in\n' ...
                       'M_{ext} = %.0f kip-in\n' ...
                       'f_{top} = %.3f ksi\n' ...
                       'f_{bot} = %.3f ksi\n' ...
                       '\x03c6 = %.2e 1/in'], ...
                       P, e, f_axial, P*e, M, f_total_top, f_total_bot, curvature);
    
    if mean(f_total_arr) < 0
        text_x = 0.98;
        text_ha = 'right';
    else
        text_x = 0.02;
        text_ha = 'left';
    end
    
    text(text_x, 0.98, info_str, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', text_ha, ...
        'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
        'Margin', 3);
    
    % Set axis limits with some padding
    curr_xlim = xlim;
    xlim([curr_xlim(1)*1.1 - 0.5, curr_xlim(2)*1.1 + 0.5]);
end

%% Plot strain distribution
if strcmp(options.plot_type, 'both') || strcmp(options.plot_type, 'strain')
    if strcmp(options.plot_type, 'both')
        subplot(1, 2, 2);
    end
    
    hold on;
    grid on;
    box on;
    
    % Convert to microstrain for readability
    epsilon_axial_mu     = epsilon_axial * 1e6;
    epsilon_prebend_mu   = epsilon_prebend * 1e6;
    epsilon_prestress_mu = epsilon_prestress * 1e6;
    epsilon_external_mu  = epsilon_external * 1e6;
    epsilon_total_mu     = epsilon_total * 1e6;
    
    h_section = max(y_from_bottom);
    
    % Plot strain components if requested
    if options.show_components
        % Prestress axial (uniform) - dashed green
        h1 = plot(epsilon_axial_mu, y_from_bottom, 'g--', 'LineWidth', 1.5, ...
            'DisplayName', 'Prestress Axial');
        
        % Prestress bending (linear) - dash-dot green
        h2 = plot(epsilon_prebend_mu, y_from_bottom, 'g-.', 'LineWidth', 1.5, ...
            'DisplayName', 'Prestress Bending');
        
        % Combined prestress - solid green
        h3 = plot(epsilon_prestress_mu, y_from_bottom, 'g-', 'LineWidth', 2, ...
            'DisplayName', 'Prestress Combined');
        
        % External load - blue
        h4 = plot(epsilon_external_mu, y_from_bottom, 'b-', 'LineWidth', 1.5, ...
            'DisplayName', 'External Load');
    end
    
    % Plot total strain - red
    h5 = plot(epsilon_total_mu, y_from_bottom, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Total');
    
    % Fill tension and compression regions
    eps_pos = max(epsilon_total_mu, 0);
    eps_neg = min(epsilon_total_mu, 0);
    
    fill([eps_pos, zeros(size(eps_pos))], [y_from_bottom, fliplr(y_from_bottom)], ...
        [1, 0.7, 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([eps_neg, zeros(size(eps_neg))], [y_from_bottom, fliplr(y_from_bottom)], ...
        [0.7, 0.7, 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Zero line
    plot([0, 0], [0, h_section], 'k-', 'LineWidth', 1);
    
    % Centroid line
    plot(xlim, [yb, yb], 'k--', 'LineWidth', 0.8);
    text(max(xlim)*0.95, yb, ' C.G.', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'right', 'FontSize', 9);
    
    % Neutral axis
    if y_na_within_section
        plot(xlim, [y_na_from_bottom, y_na_from_bottom], 'm-', 'LineWidth', 2);
        text(max(xlim)*0.95, y_na_from_bottom, sprintf(' N.A.'), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
            'FontSize', 9, 'Color', 'm', 'FontWeight', 'bold');
    end
    
    % Annotations for extreme fiber strains
    eps_top = epsilon_total_mu(end);
    eps_bot = epsilon_total_mu(1);
    
    text(eps_top, h_section, sprintf('  %.0f \x03bc\x03b5', eps_top), ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');
    text(eps_bot, 0, sprintf('  %.0f \x03bc\x03b5', eps_bot), ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');
    
    % Crushing strain limit
    eps_crush = -3000;
    plot([eps_crush, eps_crush], [0, h_section], 'k:', 'LineWidth', 1.5);
    text(eps_crush, h_section*0.02, ' \epsilon_{cu} = -3000 \mu\epsilon', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8, 'Rotation', 90);
    
    % Labels
    xlabel('Strain, \epsilon (microstrain)', 'FontSize', 11);
    ylabel('Distance from bottom (in)', 'FontSize', 11);
    title('Strain Distribution', 'FontWeight', 'bold', 'FontSize', 12);
    
    % Legend
    if options.show_components
        legend([h1, h2, h3, h4, h5], ...
            {'Prestress Axial', 'Prestress Bending', ...
             'Prestress Combined', 'External Load', 'Total'}, ...
            'Location', 'best', 'FontSize', 8);
    else
        legend(h5, {'Total'}, 'Location', 'best', 'FontSize', 9);
    end
    
    % Add text box
    info_str = sprintf(['\x03b5_{top} = %.0f \x03bc\x03b5\n' ...
                       '\x03b5_{bot} = %.0f \x03bc\x03b5\n' ...
                       'E_c = %.0f ksi'], ...
                       eps_top, eps_bot, Ec);
    
    if mean(epsilon_total_mu) < 0
        text_x = 0.98;
        text_ha = 'right';
    else
        text_x = 0.02;
        text_ha = 'left';
    end
    
    text(text_x, 0.98, info_str, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', text_ha, ...
        'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
        'Margin', 3);
    
    % Set axis limits
    curr_xlim = xlim;
    xlim([min(curr_xlim(1)*1.1, eps_crush*1.1), curr_xlim(2)*1.1 + 50]);
end

%% Add super title if both plots
if strcmp(options.plot_type, 'both')
    sgtitle(sprintf('Stress and Strain Distribution at x = %.1f in (%.1f%% of span)', ...
        x_actual, x_actual/L*100), 'FontSize', 14, 'FontWeight', 'bold');
end

%% Print summary to command window
fprintf('\n========================================\n');
fprintf('  STRESS-STRAIN DISTRIBUTION SUMMARY\n');
fprintf('========================================\n');
fprintf('Location: x = %.1f in (%.1f%% of span)\n', x_actual, x_actual/L*100);
fprintf('\nInternal Forces:\n');
fprintf('  Prestress P = %.1f kips\n', P);
fprintf('  Eccentricity e = %.2f in\n', e);
fprintf('  Prestress Moment M_p = P*e = %.0f kip-in\n', P*e);
fprintf('  External Moment M_ext = %.0f kip-in\n', M);
fprintf('  Net Moment M_net = P*e - M = %.0f kip-in\n', M_net);
fprintf('\nPrestress Stress Components:\n');
fprintf('  Axial (-P/A):                  %.4f ksi (uniform)\n', f_axial);
fprintf('  Bending at top (+P*e*yt/I):    %.4f ksi\n', f_prestress_bend_top);
fprintf('  Bending at bot (-P*e*yb/I):    %.4f ksi\n', f_prestress_bend_bot);
fprintf('  Combined prestress at top:     %.4f ksi\n', f_prestress_top);
fprintf('  Combined prestress at bot:     %.4f ksi\n', f_prestress_bot);
fprintf('\nExternal Load Stresses (-M*y/I):\n');
fprintf('  Top:    %.4f ksi\n', f_external_top);
fprintf('  Bottom: %.4f ksi\n', f_external_bot);
fprintf('\nTotal Stresses:\n');
fprintf('  Top fiber:    f = %.4f ksi (%.0f microstrain)\n', f_total_top, f_total_top/Ec*1e6);
fprintf('  Bottom fiber: f = %.4f ksi (%.0f microstrain)\n', f_total_bot, f_total_bot/Ec*1e6);
fprintf('\nNeutral Axis (where f = 0):\n');
fprintf('  y_NA from centroid = %.2f in\n', y_na_from_centroid);
fprintf('  y_NA from bottom = %.2f in\n', y_na_from_bottom);
if y_na_within_section
    fprintf('  Status: Within section\n');
else
    if isnan(y_na_from_bottom)
        fprintf('  Status: At infinity (pure axial, no net moment)\n');
    elseif y_na_from_bottom < 0
        fprintf('  Status: Below section (entire section in %s)\n', ...
            ternary(f_total_bot < 0, 'compression', 'tension'));
    else
        fprintf('  Status: Above section (entire section in %s)\n', ...
            ternary(f_total_top < 0, 'compression', 'tension'));
    end
end
fprintf('\nCurvature:\n');
fprintf('  phi = %.4e 1/in = %.4f 1/ft\n', curvature, curvature*12);
fprintf('========================================\n\n');

end

%% Helper function
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
