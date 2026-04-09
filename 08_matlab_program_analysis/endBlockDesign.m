function [eb] = endBlockDesign(beam, section, materials, prestress)
% ENDBLOCKDESIGN - End block (anchorage zone) design using Gergely-Sozen method
%
% Computes bursting and spalling reinforcement for the end block of a
% prestressed concrete beam using the linear elastic free-body approach
% (Naaman Sec. 4.17, Gergely & Sozen 1967).
%
% METHOD:
%   1. Compute elastic stresses f(y) = F_i/A + F_i*e_0*(yc-y)/Ix
%   2. Integrate f(y)*b(y)*(y_cut - y) dy to get M_concrete at each y_cut
%   3. Subtract M_prestress = F_i*(y_cut - y_ps) when y_cut > y_ps
%   4. Find M_max -> T_burst = M_max / lever_arm
%   5. Size stirrups at f_s = 20 ksi (WSD)
%
%   Two lever arm options (set 'method' variable in Section 6):
%     'deep_beam'     : lever_arm = 0.8*h  (Prof. Fafitis, CEE 530)
%     'gergely_sozen' : lever_arm = 3h/4   (Naaman Sec. 4.17)
%
% SIGN CONVENTION: compression = positive, tension = negative
% COORDINATE ORIGIN: y = 0 at bottom of stems
%
% Inputs:
%   beam      - Beam geometry (needs .L)
%   section   - Cross-section (.vertices, .A, .Ix, .yc, .yb, .yt)
%   materials - Material properties (.fci, .fc, .fy)
%   prestress - Tendon data (.tendons{i}, .losses)
%
% Output:
%   eb - Structure with all end block design results
%
% Reference:
%   Naaman, "Prestressed Concrete", 2nd ed., Sec. 4.17, Example 4.17.3
%   Gergely & Sozen, PCI Journal, 12(2): 63-75, 1967
%   ACI 318-19, Section 25.9

%% ================================================================
%%  1. EXTRACT SECTION PROPERTIES
%% ================================================================
A   = section.A;
Ix  = section.Ix;
yc  = section.yc;
yb  = section.yb;
yt  = section.yt;
h   = max(section.vertices(:,2)) - min(section.vertices(:,2));
y_min = min(section.vertices(:,2));
St  = Ix / yt;   % top section modulus
Sb  = Ix / yb;   % bottom section modulus

%% ================================================================
%%  2. COMPUTE PRESTRESS FORCE AND TENDON LOCATION AT SUPPORT
%% ================================================================
F_i   = 0;      % total initial prestress force (kip)
sum_Fy = 0;     % for force-weighted y centroid

for i = 1:length(prestress.tendons)
    t = prestress.tendons{i};
    Pi = t.Aps * t.fpi;   % initial force for this tendon (kip)
    F_i = F_i + Pi;

    % Get tendon y at support (x = 0)
    if isfield(t, 'y') && ~isempty(t.y)
        y_t = t.y(1);   % first point of computed profile
    elseif isfield(t, 'y_position') && ~isempty(t.y_position)
        y_t = t.y_position;
    elseif isfield(t, 'y_profile') && ~isempty(t.y_profile)
        y_t = t.y_profile(1);
    else
        error('Cannot determine tendon %d y-position at support.', i);
    end

    sum_Fy = sum_Fy + Pi * y_t;
end

y_ps = sum_Fy / F_i;       % resultant tendon y at support (in)
e_0  = yc - y_ps;           % eccentricity (positive when tendon below centroid)

%% ================================================================
%%  3. BUILD FINE Y-GRID WITH REFINEMENT
%% ================================================================
% Identify key y-levels for refinement
y_vertices = unique(section.vertices(:,2));
key_y = unique([y_ps, yc, y_vertices']);

% Base grid: 500 points
N_base = 500;
y_grid = linspace(y_min, y_min + h, N_base);

% Add refined points near key locations (within +/- 0.5 in, 20 extra points each)
for k = 1:length(key_y)
    yk = key_y(k);
    if yk > y_min && yk < y_min + h
        y_local = linspace(max(y_min, yk - 0.5), min(y_min + h, yk + 0.5), 20);
        y_grid = [y_grid, y_local]; %#ok<AGROW>
    end
end

% Add exact key points
y_grid = unique([y_grid, key_y]);
y_grid = sort(y_grid);
N = length(y_grid);

%% ================================================================
%%  4. COMPUTE ELASTIC STRESS f(y) AND SECTION WIDTH b(y)
%% ================================================================
% Stress from prestress only (compression = positive):
%   f(y) = F_i/A  +  M_ps * (yc - y) / Ix
% where M_ps = F_i * e_0 (hogging moment when tendon is below centroid)
% This gives:
%   At bottom (y=0): f = F_i/A + F_i*e_0*yc/Ix   (more compression)
%   At top (y=h):    f = F_i/A - F_i*e_0*yt/Ix    (less compression)
%
% Equivalently: f(y) = F_i/A + F_i*e_0*(yc - y)/Ix
f_grid = F_i/A + F_i * e_0 * (yc - y_grid) / Ix;   % stress (ksi)
b_grid = zeros(1, N);

for i = 1:N
    b_grid(i) = computeWidthAtY(section.vertices, y_grid(i));
end

% Force per unit height (kip/in)
q_grid = f_grid .* b_grid;

% Fiber stresses for reporting
f_top = F_i/A - F_i*e_0*yt/Ix;
f_bot = F_i/A + F_i*e_0*yb/Ix;

%% ================================================================
%%  5. COMPUTE NET MOMENT ON HORIZONTAL PLANES
%% ================================================================
% M_concrete(y_cut) = integral_0^y_cut of q(y) * (y_cut - y) dy
%                   = y_cut * Q(y_cut) - S(y_cut)
% where Q = integral of q dy, S = integral of q*y dy

Q = cumtrapz(y_grid, q_grid);        % cumulative force (kip)
S = cumtrapz(y_grid, q_grid .* y_grid);  % cumulative first moment (kip-in)

M_concrete = y_grid .* Q - S;         % moment from concrete stresses (kip-in)

% M_prestress: moment from concentrated prestress force at y_ps
M_prestress = zeros(1, N);
for i = 1:N
    if y_grid(i) > y_ps
        M_prestress(i) = -F_i * (y_grid(i) - y_ps);
    end
end

% Net moment
M_net = M_concrete + M_prestress;

% Find maximum
[M_max, idx_max] = max(abs(M_net));
y_Mmax = y_grid(idx_max);

%% ================================================================
%%  6. BURSTING FORCE
%% ================================================================
% Two methods available:
%   'deep_beam'       — T = M_max / (0.8*h)   (Prof. Fafitis, CEE 530)
%   'gergely_sozen'   — T = M_max / (3*h/4)   (Naaman Sec. 4.17)
%
% Deep beam: z = 0.2h from end face, lever arm = h - z = 0.8h
% Gergely-Sozen: T at h/4, C at h, lever arm = 3h/4

method = 'deep_beam';   % <-- CHANGE HERE to switch method

switch method
    case 'deep_beam'
        lever_arm = 0.8 * h;
        method_label = 'Deep Beam (Fafitis CEE 530)';
    case 'gergely_sozen'
        lever_arm = 3 * h / 4;
        method_label = 'Gergely-Sozen (Naaman Sec. 4.17)';
    otherwise
        error('Unknown method: %s. Use ''deep_beam'' or ''gergely_sozen''.', method);
end

T_burst = M_max / lever_arm;

%% ================================================================
%%  7. STIRRUP DESIGN (WSD: f_s = 20 ksi)
%% ================================================================
f_s = 20.0;   % ksi (working stress design, per Naaman / Fafitis)

% Default #3 U-stirrups
bar_size = '#3';
A_bar = 0.11;           % in^2 per leg
Av_bar = 2 * A_bar;     % in^2 per stirrup (2 legs)

As_req = T_burst / f_s;
n_stirrups = ceil(As_req / Av_bar);
spacing = h / n_stirrups;

%% ================================================================
%%  8. SPALLING REINFORCEMENT (ACI 318 Sec. 25.9.4.4.5)
%% ================================================================
P_pu = 1.2 * F_i;                        % factored prestress force (kip)
T_spalling = 0.02 * P_pu;                % spalling force (kip)
fy_s = materials.fy;                      % mild steel yield (ksi)
As_spalling = T_spalling / fy_s;          % required area (in^2)

%% ================================================================
%%  9. ACI APPROXIMATE METHOD COMPARISON
%% ================================================================
% Estimate h_anc from tendon spread at end
y_tendons_end = zeros(1, length(prestress.tendons));
for i = 1:length(prestress.tendons)
    t = prestress.tendons{i};
    if isfield(t, 'y') && ~isempty(t.y)
        y_tendons_end(i) = t.y(1);
    elseif isfield(t, 'y_position') && ~isempty(t.y_position)
        y_tendons_end(i) = t.y_position;
    elseif isfield(t, 'y_profile') && ~isempty(t.y_profile)
        y_tendons_end(i) = t.y_profile(1);
    end
end
h_anc = max(y_tendons_end) - min(y_tendons_end);
if h_anc < 1
    h_anc = 1;   % minimum 1 in if all tendons at same height
end
e_anc = abs(yc - y_ps);

T_burst_ACI = 0.25 * P_pu * (1 - h_anc / h);
d_burst_ACI = 0.5 * (h - e_anc);
phi_v = 0.75;
As_ACI = T_burst_ACI / (phi_v * fy_s);

%% ================================================================
%%  10. PACK OUTPUT STRUCT
%% ================================================================
eb.F_i         = F_i;
eb.y_ps        = y_ps;
eb.e_0         = e_0;
eb.h           = h;
eb.f_top       = f_top;
eb.f_bot       = f_bot;
eb.y_net       = y_grid;
eb.M_concrete  = M_concrete;
eb.M_prestress = M_prestress;
eb.M_net       = M_net;
eb.M_max       = M_max;
eb.y_Mmax      = y_Mmax;
eb.method      = method;
eb.method_label = method_label;
eb.lever_arm   = lever_arm;
eb.T_burst     = T_burst;
eb.T_burst_GS  = M_max / (3*h/4);    % Gergely-Sozen (for comparison)
eb.T_burst_DB  = M_max / (0.8*h);    % Deep Beam (for comparison)
eb.f_s         = f_s;
eb.As_req      = As_req;
eb.bar_size    = bar_size;
eb.Av_bar      = Av_bar;
eb.n_stirrups  = n_stirrups;
eb.spacing     = spacing;
eb.T_spalling  = T_spalling;
eb.As_spalling = As_spalling;
eb.T_burst_ACI = T_burst_ACI;
eb.d_burst_ACI = d_burst_ACI;
eb.As_ACI      = As_ACI;
eb.f_grid      = f_grid;
eb.b_grid      = b_grid;

%% ================================================================
%%  11. CONSOLE OUTPUT
%% ================================================================
printEndBlockResults(eb, section, materials);

%% ================================================================
%%  12. GENERATE FIGURE
%% ================================================================
plotEndBlockResults(eb, section);

end


%% =====================================================================
%%  LOCAL FUNCTIONS
%% =====================================================================

function b = computeWidthAtY(vertices, y)
% COMPUTEWIDTHATY - Compute total section width at height y by polygon edge intersection.
%
% For each edge of the polygon, check if the horizontal line at y intersects it.
% Collect all x-intercepts, sort them, and sum pairwise widths.

n = size(vertices, 1);
eps_y = 1e-6;   % small offset to avoid exact vertex hits

x_intercepts = [];

for i = 1:n
    j = mod(i, n) + 1;   % next vertex (wraps around)

    y1 = vertices(i, 2);
    y2 = vertices(j, 2);
    x1 = vertices(i, 1);
    x2 = vertices(j, 1);

    % Skip horizontal edges
    if abs(y2 - y1) < eps_y
        continue;
    end

    % Check if y is within the edge's y-range (exclusive of endpoints to avoid double-counting)
    y_lo = min(y1, y2);
    y_hi = max(y1, y2);

    if y >= y_lo - eps_y && y <= y_hi + eps_y
        % Linear interpolation for x at this y
        t = (y - y1) / (y2 - y1);
        if t >= -eps_y && t <= 1 + eps_y
            x_int = x1 + t * (x2 - x1);
            x_intercepts(end+1) = x_int; %#ok<AGROW>
        end
    end
end

% Sort intercepts and compute pairwise widths
x_intercepts = sort(x_intercepts);

% Remove near-duplicate intercepts (from vertex double-hits)
if length(x_intercepts) > 1
    x_clean = x_intercepts(1);
    for i = 2:length(x_intercepts)
        if abs(x_intercepts(i) - x_clean(end)) > 0.01
            x_clean(end+1) = x_intercepts(i); %#ok<AGROW>
        end
    end
    x_intercepts = x_clean;
end

% Sum pairwise widths
b = 0;
for i = 1:2:length(x_intercepts)-1
    b = b + abs(x_intercepts(i+1) - x_intercepts(i));
end

end


function printEndBlockResults(eb, section, materials)
% Print detailed end block design results to console.

fprintf('\n======================================================\n');
fprintf('   END BLOCK DESIGN -- %s\n', eb.method_label);
fprintf('   CEE 530 Prestressed Concrete (Prof. Fafitis)\n');
fprintf('======================================================\n');

fprintf('\n------ SECTION PROPERTIES ------\n');
fprintf('  h   = %.2f in\n', eb.h);
fprintf('  A   = %.3f in^2\n', section.A);
fprintf('  Ix  = %.1f in^4\n', section.Ix);
fprintf('  yc  = %.4f in (from bottom)\n', section.yc);
fprintf('  St  = %.1f in^3\n', section.Ix / section.yt);
fprintf('  Sb  = %.1f in^3\n', section.Ix / section.yb);

fprintf('\n------ PRESTRESS AT SUPPORT (x = 0) ------\n');
fprintf('  F_i  = %.1f kip (total initial, before losses)\n', eb.F_i);
fprintf('  y_ps = %.2f in (resultant tendon location from bottom)\n', eb.y_ps);
fprintf('  e_0  = %.2f in (eccentricity below centroid)\n', eb.e_0);

fprintf('\n------ STRESS DISTRIBUTION AT END FACE ------\n');
fprintf('  f(y) = F_i/A + F_i*e_0*(y - yc)/Ix\n');
fprintf('  f_top = %+.4f ksi\n', eb.f_top);
fprintf('  f_bot = %+.4f ksi\n', eb.f_bot);

fprintf('\n------ NET MOMENT ON HORIZONTAL PLANES ------\n');
fprintf('  %8s | %14s | %14s | %14s\n', 'y (in)', 'M_conc (k-in)', 'M_ps (k-in)', 'M_net (k-in)');
fprintf('  %8s-|-%14s-|-%14s-|-%14s\n', '--------', '--------------', '--------------', '--------------');

% Print at selected y-locations
y_print = unique(sort([0, 4, 8, 12, 16, eb.y_ps, 20, 24, section.yc, ...
    min(section.vertices(:,2)) + eb.h * [0.25 0.5 0.75 1.0]]));
% Keep only values within range
y_print = y_print(y_print >= min(eb.y_net) & y_print <= max(eb.y_net));

for k = 1:length(y_print)
    [~, idx] = min(abs(eb.y_net - y_print(k)));
    fprintf('  %8.2f | %14.2f | %14.2f | %14.2f', ...
        eb.y_net(idx), eb.M_concrete(idx), eb.M_prestress(idx), eb.M_net(idx));
    if abs(eb.y_net(idx) - eb.y_ps) < 0.1
        fprintf('  <-- tendon');
    end
    if abs(eb.y_net(idx) - section.yc) < 0.1
        fprintf('  <-- CGC');
    end
    if abs(eb.y_net(idx) - eb.y_Mmax) < 0.1
        fprintf('  <-- M_max');
    end
    fprintf('\n');
end

fprintf('\n------ MAXIMUM NET MOMENT ------\n');
fprintf('  M_max = %.2f kip-in at y = %.2f in\n', eb.M_max, eb.y_Mmax);

fprintf('\n------ BURSTING FORCE [%s] ------\n', eb.method_label);
fprintf('  Lever arm = %.2f * h = %.2f * %.2f = %.2f in\n', ...
    eb.lever_arm/eb.h, eb.lever_arm/eb.h, eb.h, eb.lever_arm);
fprintf('  T_burst = M_max / lever_arm = %.2f / %.2f = %.2f kip\n', ...
    eb.M_max, eb.lever_arm, eb.T_burst);

fprintf('\n  --- Method Comparison ---\n');
fprintf('  Deep Beam (Fafitis):    T = M_max / (0.8h)  = %.2f / %.2f = %.2f kip\n', ...
    eb.M_max, 0.8*eb.h, eb.T_burst_DB);
fprintf('  Gergely-Sozen (Naaman): T = M_max / (3h/4)  = %.2f / %.2f = %.2f kip\n', ...
    eb.M_max, 3*eb.h/4, eb.T_burst_GS);

fprintf('\n------ STIRRUP DESIGN (WSD: f_s = %.0f ksi) ------\n', eb.f_s);
fprintf('  A_s_req = T / f_s = %.2f / %.0f = %.4f in^2\n', eb.T_burst, eb.f_s, eb.As_req);
fprintf('  Using %s U-stirrups (A_v = %.2f in^2 per stirrup, 2 legs)\n', eb.bar_size, eb.Av_bar);
fprintf('  n = ceil(%.4f / %.2f) = %d stirrups\n', eb.As_req, eb.Av_bar, eb.n_stirrups);
fprintf('  Spacing = h / n = %.2f / %d = %.1f in\n', eb.h, eb.n_stirrups, eb.spacing);
fprintf('  Zone: 0 to %.0f in from end face\n', eb.h);

fprintf('\n------ SPALLING REINFORCEMENT (ACI 318 Sec. 25.9.4.4.5) ------\n');
fprintf('  P_pu = 1.2 * F_i = 1.2 * %.1f = %.1f kip\n', eb.F_i, 1.2*eb.F_i);
fprintf('  T_spalling = 0.02 * P_pu = 0.02 * %.1f = %.2f kip\n', 1.2*eb.F_i, eb.T_spalling);
fprintf('  A_s_spalling = T_spalling / f_y = %.2f / %.0f = %.4f in^2\n', ...
    eb.T_spalling, materials.fy, eb.As_spalling);

fprintf('\n------ ACI 318 APPROXIMATE METHOD (comparison) ------\n');
fprintf('  T_burst_ACI = 0.25 * P_pu * (1 - h_anc/h) = %.2f kip\n', eb.T_burst_ACI);
fprintf('  d_burst_ACI = 0.5 * (h - e_anc) = %.2f in\n', eb.d_burst_ACI);
fprintf('  A_s_ACI = T_burst_ACI / (phi*f_y) = %.2f / (0.75*%.0f) = %.4f in^2\n', ...
    eb.T_burst_ACI, materials.fy, eb.As_ACI);

fprintf('\n======================================================\n');
end


function plotEndBlockResults(eb, section)
% Generate 5 separate figures for end block design:
%   1. Stress distribution f(y) vs y
%   2. Force distribution q(y) = f(y)*b(y) vs y
%   3. Net moment on horizontal planes M_net(y) vs y
%   4. Reinforcement layout (elevation)
%   (Cross-section and stress/strain at x=0 are plotted by the caller
%    using the shared plotSection / plotSectionStressStrain templates.)

y = eb.y_net;
y_min_sec = min(section.vertices(:,2));

%% Figure: Stress distribution f(y) vs y
figure('Name', 'EndBlock - Stress Distribution', ...
    'Position', [100 100 500 600]);
hold on;

% Stress line
plot(eb.f_grid, y, 'b-', 'LineWidth', 2.5);

% Shade stress area
fill([eb.f_grid, 0, 0], [y, y(end), y(1)], ...
    [0.7 0.8 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(eb.f_grid, y, 'b-', 'LineWidth', 2.5);

% Zero line
plot([0 0], [y_min_sec, y_min_sec + eb.h], 'k-', 'LineWidth', 0.8);

% Annotate fiber stresses
text(eb.f_top, y_min_sec + eb.h, sprintf('  %+.1f psi', eb.f_top*1000), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b', 'VerticalAlignment', 'middle');
text(eb.f_bot, y_min_sec, sprintf('  %+.1f psi', eb.f_bot*1000), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b', 'VerticalAlignment', 'middle');

% Mark tendon and CGC
yline(eb.y_ps, ':', sprintf('y_{ps} = %.1f in', eb.y_ps), ...
    'Color', 'r', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
yline(section.yc, ':', sprintf('y_c = %.1f in', section.yc), ...
    'Color', [0.3 0.3 0.3], 'LineWidth', 1, 'LabelHorizontalAlignment', 'left', 'FontSize', 9);

xlabel('Stress, f(y) (ksi)', 'FontSize', 12);
ylabel('Height, y (in)', 'FontSize', 12);
title('End Face Stress Distribution', 'FontSize', 14);
grid on;
hold off;

%% Figure: Force distribution q(y) = f(y)*b(y) vs y
figure('Name', 'EndBlock - Force Distribution', ...
    'Position', [620 100 500 600]);
hold on;

q_grid = eb.f_grid .* eb.b_grid;   % kip/in

% Force line
plot(q_grid, y, 'b-', 'LineWidth', 2.5);

% Shade force area
fill([q_grid, 0, 0], [y, y(end), y(1)], ...
    [0.7 0.8 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(q_grid, y, 'b-', 'LineWidth', 2.5);

% Zero line
plot([0 0], [y_min_sec, y_min_sec + eb.h], 'k-', 'LineWidth', 0.8);

% Annotate key values at section transitions
y_vertices = unique(sort(section.vertices(:,2)));
y_label = unique(sort([y_min_sec, y_vertices', y_min_sec + eb.h, eb.y_ps, section.yc]));
y_label = y_label(y_label >= min(y) & y_label <= max(y));
for k = 1:length(y_label)
    [~, idx_k] = min(abs(y - y_label(k)));
    q_val = q_grid(idx_k);
    if abs(q_val) > max(abs(q_grid)) * 0.02   % skip near-zero
        plot(q_val, y(idx_k), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        text(q_val, y(idx_k), sprintf('  %.0f lb/in (y=%.0f)', q_val*1000, y(idx_k)), ...
            'FontSize', 8, 'VerticalAlignment', 'middle');
    end
end

% Mark tendon and CGC
yline(eb.y_ps, ':', 'Color', 'r', 'LineWidth', 1);
yline(section.yc, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);

xlabel('Force per unit height, q(y) = f(y) \cdot b(y) (kip/in)', 'FontSize', 11);
ylabel('Height, y (in)', 'FontSize', 12);
title('Force Distribution Along Section Depth', 'FontSize', 14);
grid on;
hold off;

%% Figure: Net moment diagram
figure('Name', 'EndBlock - Net Moment Diagram', ...
    'Position', [720 100 600 500]);
hold on;

% Plot all three moment curves
plot(eb.M_concrete, y, 'b--', 'LineWidth', 1.5);
plot(eb.M_prestress, y, 'r--', 'LineWidth', 1.5);
plot(eb.M_net, y, 'k-', 'LineWidth', 2.5);

% Mark M_max
plot(eb.M_max * sign(eb.M_net(abs(eb.M_net) == eb.M_max)), eb.y_Mmax, ...
    'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

% Mark tendon and CGC
yline(eb.y_ps, ':', 'y_{ps}', 'Color', 'r', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
yline(section.yc, ':', 'y_c', 'Color', [0.3 0.3 0.3], 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');

xlabel('Moment (kip-in)', 'FontSize', 12);
ylabel('Height, y (in)', 'FontSize', 12);
title('Net Moment on Horizontal Planes', 'FontSize', 14);
legend('M_{concrete}', 'M_{prestress}', 'M_{net}', ...
    sprintf('M_{max} = %.0f kip-in at y = %.1f in', eb.M_max, eb.y_Mmax), ...
    'Location', 'best', 'FontSize', 10);
grid on;
hold off;

%% Figure: Reinforcement layout (elevation view of end block)
figure('Name', 'EndBlock - Reinforcement Layout', ...
    'Position', [100 100 700 500]);
hold on;

% Draw section outline (elevation: x = distance from end, y = height)
rect_x = [0, eb.h, eb.h, 0, 0];
rect_y = [y_min_sec, y_min_sec, y_min_sec + eb.h, y_min_sec + eb.h, y_min_sec];
plot(rect_x, rect_y, 'k-', 'LineWidth', 1.5);

% Draw tendon
plot([0, eb.h], [eb.y_ps, eb.y_ps], 'r-', 'LineWidth', 2);
plot(0, eb.y_ps, 'r>', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(eb.h * 0.02, eb.y_ps + 1, sprintf('F_i = %.0f k', eb.F_i), ...
    'Color', 'r', 'FontSize', 10);

% Draw stirrups
stirrup_x = linspace(eb.spacing/2, eb.h - eb.spacing/2, eb.n_stirrups);
for k = 1:length(stirrup_x)
    plot([stirrup_x(k), stirrup_x(k)], [y_min_sec + 1.5, y_min_sec + eb.h - 1.5], ...
        'b-', 'LineWidth', 2);
end

% Draw CGC
plot([0, eb.h], [section.yc, section.yc], 'k--', 'LineWidth', 0.8);
text(eb.h * 0.6, section.yc + 1, 'CGC', 'FontSize', 10);

% Draw M_max location
yline(eb.y_Mmax, '-.', 'Color', [0.8 0 0], 'LineWidth', 1);
text(eb.h * 0.6, eb.y_Mmax - 1.5, sprintf('M_{max} at y = %.1f in', eb.y_Mmax), ...
    'Color', [0.8 0 0], 'FontSize', 9);

title(sprintf('End Block Reinforcement Layout\n%d %s U-stirrups @ %.1f in (zone = %.0f in from end)', ...
    eb.n_stirrups, eb.bar_size, eb.spacing, eb.h), 'FontSize', 14);
xlabel('Distance from end face (in)', 'FontSize', 12);
ylabel('Height, y (in)', 'FontSize', 12);
axis equal;
xlim([-2, eb.h + 2]);
ylim([y_min_sec - 2, y_min_sec + eb.h + 4]);
grid on;
hold off;

end
