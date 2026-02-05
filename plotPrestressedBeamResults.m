function plotPrestressedBeamResults(results, options)
% PLOTPRESTRESSEDBEAMRESULTS - Plot analysis results for prestressed beam
%
% Sign Convention (Universal Beam Convention):
%   +M = Sagging (tension at bottom)
%   -M = Hogging (tension at top)
%   M_external: gravity loads -> positive (sagging)
%   M_prestress: tendon below CG -> negative (hogging/camber)
%   M_net = M_external + M_prestress
%
% Inputs:
%   results - Analysis results structure from analyzePrestressedBeam
%   options - Optional structure for plot customization

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

%% ========================================================================
%  BEAM SCHEMATIC - Bottom surface at y = 0
%% ========================================================================
function plotBeamSchematic(results)

x = results.x;
L = results.L;
section = results.section;
prestress = results.prestress;

hold on;
grid on;
box on;

% Beam dimensions - bottom at y = 0
h = section.yt + section.yb;
yc = section.yc;
y_bot = 0;
y_top = h;

% Draw beam outline
plot([0, L], [yc, yc], 'k--', 'LineWidth', 1);       % Centroid (dashed)
plot([0, L], [y_bot, y_bot], 'k-', 'LineWidth', 1);   % Bottom
plot([0, L], [y_top, y_top], 'k-', 'LineWidth', 1);   % Top
plot([0, 0], [y_bot, y_top], 'k-', 'LineWidth', 1);   % Left end
plot([L, L], [y_bot, y_top], 'k-', 'LineWidth', 1);   % Right end

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
    y_tendon = yc - tendon.e;   % Actual y from bottom

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

% Supports at y = 0
support_h = 0.08 * h;
if isfield(results.loads, 'support_type')
    switch results.loads.support_type
        case 'simple'
            drawTriangle(results.reactions.x_left, y_bot, support_h, 'pin');
            drawTriangle(results.reactions.x_right, y_bot, support_h, 'roller');
        case 'cantilever'
            drawFixedSupport(results.reactions.x_support, y_bot, support_h);
    end
end

% Legend for tendons
legend_entries = cell(1, length(prestress.tendons));
for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};
    legend_entries{i} = sprintf('Tendon %d: A_{ps}=%.2f in^2, %s', ...
        i, tendon.Aps, tendon.bonding.type);
end
legend(legend_entries, 'Location', 'northeast');

xlabel('Distance along beam (in)');
ylabel('Elevation from bottom (in)');
xlim([-0.05*L, 1.1*L]);
ylim([y_bot - 0.2*h, y_top + 0.1*h]);

end

%% ========================================================================
%  AXIAL FORCE DIAGRAM
%% ========================================================================
function plotAxialForceDiagram(results, options)

x = results.x;
N = results.N;

hold on;
grid on;
box on;

fill([x, fliplr(x)], [N, zeros(size(N))], [0.2, 0.6, 1], ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(x, N, 'b-', 'LineWidth', 1.5);
plot([x(1), x(end)], [0, 0], 'k-', 'LineWidth', 0.5);

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

text(0.02, 0.98, 'Compression (-)', 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.5, 0.5, 0.5]);

end

%% ========================================================================
%  SHEAR FORCE DIAGRAM
%% ========================================================================
function plotShearForceDiagram(results, options)

x = results.x;
V = results.V;

hold on;
grid on;
box on;

V_pos = max(V, 0);
V_neg = min(V, 0);

fill([x, fliplr(x)], [V_pos, zeros(size(V_pos))], [1, 0.4, 0.4], ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');
fill([x, fliplr(x)], [V_neg, zeros(size(V_neg))], [0.4, 0.4, 1], ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');

plot(x, V, 'k-', 'LineWidth', 1.5);
plot([x(1), x(end)], [0, 0], 'k-', 'LineWidth', 0.5);

[V_max, idx_max] = max(V);
[V_min, idx_min] = min(V);

text(x(idx_max), V_max, sprintf('%.1f kips', V_max), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(x(idx_min), V_min, sprintf('%.1f kips', V_min), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

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

%% ========================================================================
%  BENDING MOMENT DIAGRAM - CORRECTED SIGN CONVENTION
%
%  M_external:  from gravity loads -> positive (sagging)
%  M_prestress: from P*e with tendon below CG -> negative (hogging)
%  M_net = M_external + M_prestress
%
%  Y-axis: positive up (sagging up, hogging down)
%% ========================================================================
function plotMomentDiagram(results, options)

x = results.x;
M = results.M;                      % External moment: +M = sagging
M_prestress = results.M_prestress;  % Prestress moment: -M = hogging (tendon below CG)

hold on;
grid on;
box on;

% Net moment
M_net = M + M_prestress;

% Color fill for external moment
M_sag = max(M, 0);
M_hog = min(M, 0);

fill([x, fliplr(x)], [M_sag, zeros(size(M_sag))], [1, 0.6, 0.6], ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([x, fliplr(x)], [M_hog, zeros(size(M_hog))], [0.6, 0.6, 1], ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot lines
h1 = plot(x, M, 'k-', 'LineWidth', 2);
h2 = plot(x, M_prestress, 'g--', 'LineWidth', 1.5);
h3 = plot(x, M_net, 'r-', 'LineWidth', 1.5);

% Zero line
plot([x(1), x(end)], [0, 0], 'k-', 'LineWidth', 0.5);

% Annotate max external moment
[~, idx_max] = max(abs(M));
text(x(idx_max), M(idx_max), ...
    sprintf('M_{ext} = %.0f kip-in\n(%.0f kip-ft)', M(idx_max), M(idx_max)/12), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'BackgroundColor', 'w');

% Annotate prestress moment at midspan
mid_idx = round(length(x)/2);
text(x(mid_idx), M_prestress(mid_idx), ...
    sprintf('M_{ps} = %.0f kip-in', M_prestress(mid_idx)), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontSize', 8, 'Color', [0, 0.5, 0]);

% Capacity line if available
if isfield(results, 'capacity')
    phi_Mn = results.capacity.phi_Mn;
    plot([x(1), x(end)], [phi_Mn, phi_Mn], 'm--', 'LineWidth', 1);
    text(x(end), phi_Mn, sprintf('\\phi M_n = %.0f kip-in', phi_Mn), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
        'Color', 'm');
end

legend([h1, h2, h3], {'M_{external}', 'M_{prestress}', 'M_{net}'}, ...
    'Location', 'best');

xlabel('Distance (in)');
ylabel('Bending Moment, M (kip-in)');
title('Bending Moment Diagram (+M = Sagging, -M = Hogging)', 'FontWeight', 'bold');
xlim([0, results.L]);

% Sign convention note
text(0.02, 0.98, '+M = Sagging (tension at bottom)', ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Color', [0.5, 0.5, 0.5]);
text(0.02, 0.93, '-M = Hogging (tension at top)', ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Color', [0.5, 0.5, 0.5]);

end

%% ========================================================================
%  STRESS DISTRIBUTION
%% ========================================================================
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

h1 = plot(x, f_prestress, 'g-', 'LineWidth', 1.5);
h2 = plot(x, f_load, 'b-', 'LineWidth', 1.5);
h3 = plot(x, f_total, 'r-', 'LineWidth', 2);

plot([x(1), x(end)], [0, 0], 'k-', 'LineWidth', 0.5);

plot([x(1), x(end)], [stresses.fc_allow_compression, stresses.fc_allow_compression], ...
    'k--', 'LineWidth', 1);
plot([x(1), x(end)], [stresses.fc_allow_tension, stresses.fc_allow_tension], ...
    'k--', 'LineWidth', 1);

legend([h1, h2, h3], {'Prestress', 'External Load', 'Total'}, 'Location', 'best');

xlabel('Distance (in)');
ylabel('Stress (ksi)');
title(sprintf('Stress at %s', fiber_name), 'FontWeight', 'bold');
xlim([0, results.L]);

text(0.02, 0.02, 'Compression (-) / Tension (+)', 'Units', 'normalized', ...
    'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', [0.5, 0.5, 0.5]);

end

%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================
function drawTriangle(x, y, h, type)

w = h * 0.8;
vertices = [x, y; x - w/2, y - h; x + w/2, y - h];
fill(vertices(:,1), vertices(:,2), [0.3, 0.3, 0.3]);

if strcmp(type, 'roller')
    theta = linspace(0, 2*pi, 20);
    r = h * 0.15;
    for offset = [-w/4, w/4]
        cx = x + offset;
        cy = y - h - r;
        fill(cx + r*cos(theta), cy + r*sin(theta), 'w', 'EdgeColor', 'k');
    end
    plot([x - w/2, x + w/2], [y - h - 2*r, y - h - 2*r], 'k-', 'LineWidth', 2);
else
    plot([x - w/2, x + w/2], [y - h, y - h], 'k-', 'LineWidth', 2);
end

end

function drawFixedSupport(x, y, h)

w = h * 0.3;
fill([x - w, x, x, x - w], [y - h, y - h, y + h, y + h], [0.5, 0.5, 0.5]);

for i = 1:5
    yi = y - h + (i-1) * 2*h/5;
    plot([x - w, x], [yi, yi + h/5], 'k-', 'LineWidth', 0.5);
end

end
