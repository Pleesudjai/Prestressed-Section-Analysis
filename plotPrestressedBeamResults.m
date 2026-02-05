function plotPrestressedBeamResults(results, options)
% PLOTPRESTRESSEDBEAMRESULTS - Plot analysis results for prestressed beam
%
% Inputs:
%   results - Analysis results structure from analyzePrestressedBeam
%   options - Optional structure for plot customization
%
% Usage:
%   plotPrestressedBeamResults(results)
%   plotPrestressedBeamResults(results, options)

if nargin < 2
    options = struct();
end

% Default options
if ~isfield(options, 'units'), options.units = 'kip-in'; end
if ~isfield(options, 'show_prestress'), options.show_prestress = true; end
if ~isfield(options, 'show_stresses'), options.show_stresses = true; end
if ~isfield(options, 'show_section'), options.show_section = true; end
if ~isfield(options, 'colormap'), options.colormap = 'default'; end

%% Create figure with multiple subplots
fig = figure('Name', 'Prestressed Beam Analysis Results', ...
    'Position', [50, 50, 1400, 900], 'Color', 'w');

% Layout: 2 columns, variable rows
if options.show_section
    subplot_rows = 4;
else
    subplot_rows = 3;
end

%% Plot 1: Beam schematic with tendon profile
subplot(subplot_rows, 2, [1, 2]);
plotBeamSchematic(results);
title('Beam Schematic with Tendon Profile', 'FontSize', 12, 'FontWeight', 'bold');

%% Plot 2: Axial Force Diagram (N)
subplot(subplot_rows, 2, 3);
plotAxialForceDiagram(results, options);

%% Plot 3: Shear Force Diagram (V)
subplot(subplot_rows, 2, 4);
plotShearForceDiagram(results, options);

%% Plot 4: Bending Moment Diagram (M)
subplot(subplot_rows, 2, [5, 6]);
plotMomentDiagram(results, options);

%% Plot 5-6: Stress distribution (if enabled)
if options.show_stresses && subplot_rows >= 4
    subplot(subplot_rows, 2, 7);
    plotStressDistribution(results, 'top');
    
    subplot(subplot_rows, 2, 8);
    plotStressDistribution(results, 'bottom');
end

%% Adjust layout
sgtitle('Prestressed Concrete Beam Analysis', 'FontSize', 14, 'FontWeight', 'bold');

end

%% BEAM SCHEMATIC
function plotBeamSchematic(results)


x = results.x;
L = results.L;
section = results.section;
prestress = results.prestress;

hold on;
grid on;
box on;

% Beam dimensions - bottom at y=0
h = section.yt + section.yb;  % Total height
yc = section.yc;              % Centroid from bottom
y_bot = 0;
y_top = h;

% Draw beam outline
plot([0, L], [yc, yc], 'k--', 'LineWidth', 1);      % Centroid (dashed)
plot([0, L], [y_bot, y_bot], 'k-', 'LineWidth', 1);  % Bottom at y=0
plot([0, L], [y_top, y_top], 'k-', 'LineWidth', 1);  % Top
plot([0, 0], [y_bot, y_top], 'k-', 'LineWidth', 1);  % Left end
plot([L, L], [y_bot, y_top], 'k-', 'LineWidth', 1);  % Right end

% Fill beam
fill([0, L, L, 0], [y_bot, y_bot, y_top, y_top], ...
    [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Centroid label
text(L * 1.02, yc, 'C.G.', 'FontSize', 8, 'Color', [0.4, 0.4, 0.4], ...
    'VerticalAlignment', 'middle');

% Plot tendon profiles
colors = lines(length(prestress.tendons));
for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};
    y_tendon = yc - tendon.e;  % Actual y from bottom (no scale)
    
    bonded = tendon.bonded;
    diff_bonded = [true, diff(bonded) ~= 0];
    segment_starts = find(diff_bonded);
    segment_ends = [segment_starts(2:end) - 1, length(bonded)];
    
    for j = 1:length(segment_starts)
        idx = segment_starts(j):segment_ends(j);
        if bonded(segment_starts(j))
            plot(x(idx), y_tendon(idx), '-', 'Color', colors(i,:), 'LineWidth', 2);
        else
            plot(x(idx), y_tendon(idx), '--', 'Color', colors(i,:), 'LineWidth', 2);
        end
    end
end

% Supports at y=0 (bottom of beam)
support_h = 0.05 * h;
if isfield(results.loads, 'support_type')
    switch results.loads.support_type
        case 'simple'
            drawTriangle(results.reactions.x_left, y_bot, support_h, 'pin');
            drawTriangle(results.reactions.x_right, y_bot, support_h, 'roller');
        case 'cantilever'
            drawFixedSupport(results.reactions.x_support, y_bot, support_h);
    end
end

% Legend
legend_entries = cell(1, length(prestress.tendons));
for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};
    legend_entries{i} = sprintf('Tendon %d: A_{ps}=%.2f in², %s', ...
        i, tendon.Aps, tendon.bonding.type);
end
legend(legend_entries, 'Location', 'northeast');

xlabel('Distance along beam (in)');
ylabel('Elevation from bottom (in)');
xlim([-0.05*L, 1.1*L]);
ylim([y_bot - 0.15*h, y_top + 0.1*h]);

end

%% AXIAL FORCE DIAGRAM
function plotAxialForceDiagram(results, options)

x = results.x;
N = results.N;

hold on;
grid on;
box on;

% Fill area
fill([x, fliplr(x)], [N, zeros(size(N))], [0.2, 0.6, 1], ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');

% Plot line
plot(x, N, 'b-', 'LineWidth', 1.5);

% Zero line
plot([x(1), x(end)], [0, 0], 'k-', 'LineWidth', 0.5);

% Add max/min annotations
[N_max, idx_max] = max(N);
[N_min, idx_min] = min(N);

if abs(N_max) > 1e-6
    text(x(idx_max), N_max, sprintf('%.1f kips', N_max), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
if abs(N_min) > 1e-6
    text(x(idx_min), N_min, sprintf('%.1f kips', N_min), ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
end

xlabel('Distance (in)');
ylabel('Axial Force, N (kips)');
title('Axial Force Diagram', 'FontWeight', 'bold');
xlim([0, results.L]);

% Add note about sign convention
text(0.02, 0.98, 'Compression (-)', 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.5, 0.5, 0.5]);

end

%% SHEAR FORCE DIAGRAM
function plotShearForceDiagram(results, options)

x = results.x;
V = results.V;

hold on;
grid on;
box on;

% Separate positive and negative regions for coloring
V_pos = max(V, 0);
V_neg = min(V, 0);

% Fill positive shear (red)
fill([x, fliplr(x)], [V_pos, zeros(size(V_pos))], [1, 0.4, 0.4], ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');

% Fill negative shear (blue)
fill([x, fliplr(x)], [V_neg, zeros(size(V_neg))], [0.4, 0.4, 1], ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');

% Plot line
plot(x, V, 'k-', 'LineWidth', 1.5);

% Zero line
plot([x(1), x(end)], [0, 0], 'k-', 'LineWidth', 0.5);

% Add max/min annotations
[V_max, idx_max] = max(V);
[V_min, idx_min] = min(V);

text(x(idx_max), V_max, sprintf('%.1f kips', V_max), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(x(idx_min), V_min, sprintf('%.1f kips', V_min), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% If prestress shear is significant, show it
if max(abs(results.V_prestress)) > 0.1 * max(abs(V))
    plot(x, results.V_prestress, 'g--', 'LineWidth', 1);
    legend('V_{positive}', 'V_{negative}', 'V_{total}', 'V_{prestress}', ...
        'Location', 'best');
end

xlabel('Distance (in)');
ylabel('Shear Force, V (kips)');
title('Shear Force Diagram', 'FontWeight', 'bold');
xlim([0, results.L]);

end

%% BENDING MOMENT DIAGRAM
function plotMomentDiagram(results, options)

x = results.x;
M = results.M;
M_prestress = results.M_prestress;

hold on;
grid on;
box on;

% Plot with M positive downward (conventional for beams)
M_plot = -M;  % Flip for conventional plot
M_prestress_plot = M_prestress;  % Secondary moment (opposite direction)

% Separate positive and negative regions
M_pos = max(M_plot, 0);
M_neg = min(M_plot, 0);

% Fill positive moment (tension at bottom - red)
fill([x, fliplr(x)], [M_pos, zeros(size(M_pos))], [1, 0.6, 0.6], ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');

% Fill negative moment (tension at top - blue)
fill([x, fliplr(x)], [M_neg, zeros(size(M_neg))], [0.6, 0.6, 1], ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');

% Plot external moment
h1 = plot(x, M_plot, 'k-', 'LineWidth', 2);

% Plot prestress moment if enabled
h2 = plot(x, M_prestress_plot, 'g--', 'LineWidth', 1.5);

% Plot combined (net) moment
M_net = M_plot + M_prestress_plot;
h3 = plot(x, M_net, 'r-', 'LineWidth', 1.5);

% Zero line
plot([x(1), x(end)], [0, 0], 'k-', 'LineWidth', 0.5);

% Add max annotations
[M_max, idx_max] = max(abs(M));
text(x(idx_max), M_plot(idx_max), sprintf('M_{ext} = %.0f kip-in\n(%.0f kip-ft)', ...
    M(idx_max), M(idx_max)/12), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'BackgroundColor', 'w');

% Add capacity line if available
if isfield(results, 'capacity')
    phi_Mn = results.capacity.phi_Mn;
    plot([x(1), x(end)], [-phi_Mn, -phi_Mn], 'm--', 'LineWidth', 1);
    text(x(end), -phi_Mn, sprintf('\\phi M_n = %.0f kip-in', phi_Mn), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
        'Color', 'm');
end

legend([h1, h2, h3], {'M_{external}', 'M_{prestress}', 'M_{net}'}, ...
    'Location', 'best');

xlabel('Distance (in)');
ylabel('Bending Moment (kip-in)');
title('Bending Moment Diagram (Positive = Tension at Bottom)', 'FontWeight', 'bold');
xlim([0, results.L]);

% Add note
text(0.02, 0.02, 'Moment plotted on tension side', 'Units', 'normalized', ...
    'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', [0.5, 0.5, 0.5]);

end

%% STRESS DISTRIBUTION
function plotStressDistribution(results, fiber)

x = results.x;
stresses = results.stresses;

hold on;
grid on;
box on;

if strcmp(fiber, 'top')
    f_total = stresses.f_top_total;
    f_prestress = stresses.f_top_prestress;
    f_load = stresses.f_top;
    fiber_name = 'Top Fiber';
else
    f_total = stresses.f_bot_total;
    f_prestress = stresses.f_bot_prestress;
    f_load = stresses.f_bot;
    fiber_name = 'Bottom Fiber';
end

% Plot stresses
h1 = plot(x, f_prestress, 'g-', 'LineWidth', 1.5);
h2 = plot(x, f_load, 'b-', 'LineWidth', 1.5);
h3 = plot(x, f_total, 'r-', 'LineWidth', 2);

% Zero line
plot([x(1), x(end)], [0, 0], 'k-', 'LineWidth', 0.5);

% Allowable limits
plot([x(1), x(end)], [stresses.fc_allow_compression, stresses.fc_allow_compression], ...
    'k--', 'LineWidth', 1);
plot([x(1), x(end)], [stresses.fc_allow_tension, stresses.fc_allow_tension], ...
    'k--', 'LineWidth', 1);

legend([h1, h2, h3], {'Prestress', 'External Load', 'Total'}, 'Location', 'best');

xlabel('Distance (in)');
ylabel('Stress (ksi)');
title(sprintf('Stress at %s', fiber_name), 'FontWeight', 'bold');
xlim([0, results.L]);

% Add compression/tension labels
text(0.02, 0.02, 'Compression (-) / Tension (+)', 'Units', 'normalized', ...
    'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', [0.5, 0.5, 0.5]);

end

%% HELPER FUNCTIONS
function drawTriangle(x, y, h, type)
% Draw support triangle

w = h * 0.8;
vertices = [x, y; x - w/2, y - h; x + w/2, y - h];
fill(vertices(:,1), vertices(:,2), [0.3, 0.3, 0.3]);

if strcmp(type, 'roller')
    % Add circles for roller
    theta = linspace(0, 2*pi, 20);
    r = h * 0.15;
    for offset = [-w/4, w/4]
        cx = x + offset;
        cy = y - h - r;
        fill(cx + r*cos(theta), cy + r*sin(theta), 'w', 'EdgeColor', 'k');
    end
    % Ground line
    plot([x - w/2, x + w/2], [y - h - 2*r, y - h - 2*r], 'k-', 'LineWidth', 2);
else
    % Ground line for pin
    plot([x - w/2, x + w/2], [y - h, y - h], 'k-', 'LineWidth', 2);
end

end

function drawFixedSupport(x, y, h)
% Draw fixed support (wall)

w = h * 0.3;

% Wall
fill([x - w, x, x, x - w], [y - h, y - h, y + h, y + h], [0.5, 0.5, 0.5]);

% Hatching
for i = 1:5
    yi = y - h + (i-1) * 2*h/5;
    plot([x - w, x], [yi, yi + h/5], 'k-', 'LineWidth', 0.5);
end

end
