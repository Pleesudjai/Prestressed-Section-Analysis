function [beam, section, materials, prestress, reinforcement, loads] = inputPrestressedBeam_HW1()
% INPUTPRESTRESSEDBEAM_HW1 - Input data for homework 1 prestressed concrete beam
% Creates axial force, shear force, and bending moment diagrams
%
% Homework 1 Specifications:
%   - I-beam section centered at origin
%   - Single prestressing tendon at (0, 6) fully bonded
%   - P = 400 kips
%   - fc' = 6000 psi = 6 ksi
%   - Ec = 4700 ksi

%% SECTION 1: BEAM GEOMETRY
beam.L = 100;                   % Total beam length (in)
beam.num_segments = 100;        % Number of segments for analysis
beam.x = linspace(0, beam.L, beam.num_segments + 1);

%% SECTION 2: CROSS-SECTION GEOMETRY
section.vertices = [
    -8,   0;        % Bottom-left of bottom flange
    8,    0;        % Bottom-right of bottom flange
    8,    8;        % Top-right of bottom flange
    4,    8;        % Bottom-right of web
    4,    52;       % Top-right of web
    35,   52;       % Bottom-right of top flange
    35,   60;       % Top-right of top flange
    -35,  60;       % Top-left of top flange
    -35,  52;       % Bottom-left of top flange
    -4,   52;       % Top-left of web
    -4,   8;        % Bottom-left of web
    -8,   8;        % Top-left of bottom flange
];

section.A = [];
section.Ix = [];
section.Iy = [];
section.yc = [];
section.xc = [];
section.yt = [];
section.yb = [];

%% SECTION 3: MATERIAL PROPERTIES
materials.fc = 6.0;
materials.Ec = 4700;
materials.fci = 4.5;
materials.fr = 7.5 * sqrt(materials.fc * 1000) / 1000;

materials.fpu = 270;
materials.fpy = 0.9 * materials.fpu;
materials.Eps = 28500;

materials.fy = 60;
materials.Es = 29000;

%% SECTION 4: PRESTRESSING TENDONS
fpi = 0.75 * materials.fpu;
Aps_required = 400 / fpi;

% =====================================================================
% TENDON PROFILE TYPES:
%
%   'straight'         - Constant y-position along beam
%                        Requires: y_position
%
%   'linear'           - Linear eccentricity from e_start to e_end
%                        Requires: e_start, e_end
%
%   'parabolic'        - Parabolic eccentricity from 3 eccentricity values
%                        Requires: e_start, e_mid, e_end
%
%   'parabolic_y'      - Parabolic y(x) = a*x^2 + b*x + c (direct coeffs)
%                        Requires: y_coeffs = [a, b, c]
%
%   'parabolic_3pt'    - Parabolic y fitted through y_start, y_mid, y_end
%                        Requires: y_start, y_mid, y_end
%
%   'parabolic_sym'    - SYMMETRIC parabola from midspan:
%                          y(x) = a*(x - L/2)^2 + y_mid
%                        The coefficient 'a' is auto-calculated from cover:
%                          y_end = cover_bot  =>  a = (cover_bot - y_mid) / (L/2)^2
%                        Requires: y_mid, cover_bot (or y_end_override)
%                        Ensures: y(0) = y(L) = cover_bot
%
%   'harped'           - Two linear segments meeting at drape point
%                        Requires: e_start, e_mid, e_end, drape_point
%
%   'custom'           - User-defined profile vectors
%                        Requires: x_profile, e_profile
%
% CONCRETE COVER:
%   cover_bot - Clear cover from bottom fiber to tendon center (in)
%   cover_top - Clear cover from top fiber to tendon center (in)
%   These are used for:
%     (a) 'parabolic_sym' to compute end y-position automatically
%     (b) Boundary check warnings for ALL profile types
% =====================================================================

% --- Concrete cover requirements ---
tendon1.cover_bot = 3.0;        % Min cover from bottom fiber (in)
tendon1.cover_top = 3.0;        % Min cover from top fiber (in)

% --- Tendon properties ---
tendon1.Aps = Aps_required;
tendon1.fpi = fpi;

% --- ACTIVE PROFILE (uncomment ONE block) ---

% % OPTION A: Straight tendon at constant y
% tendon1.profile_type = 'straight';
% tendon1.y_position = 6;

% % OPTION B: Symmetric parabola from midspan (RECOMMENDED for draped)
%   y(x) = a*(x - L/2)^2 + y_mid
%   where a = (cover_bot - y_mid) / (L/2)^2
%   Tendon is lowest at midspan, rises to cover_bot at ends
tendon1.profile_type = 'parabolic_sym';
tendon1.y_mid = 6;              % y at midspan (lowest point, in)
  % y at ends will be = cover_bot (set above)
  % To override end position instead of using cover:
  tendon1.y_end_override = 40.80;  % optional: use this y at ends instead of cover_bot

% % OPTION C: Parabolic y from coefficients y = a*x^2 + b*x + c
% tendon1.profile_type = 'parabolic_y';
% tendon1.y_coeffs = [-0.0016, 0.16, 6];

% % OPTION D: Parabolic y from 3 control points
% tendon1.profile_type = 'parabolic_3pt';
% tendon1.y_start = 6;
% tendon1.y_mid   = 10;
% tendon1.y_end   = 6;

% % OPTION E: Parabolic eccentricity
% tendon1.profile_type = 'parabolic';
% tendon1.e_start = 30;
% tendon1.e_mid   = 35;
% tendon1.e_end   = 30;

% % OPTION F: Harped profile
% tendon1.profile_type = 'harped';
% tendon1.e_start = 10;
% tendon1.e_mid   = 35;
% tendon1.e_end   = 10;
% tendon1.drape_point = 0.5;

% Defaults for unused fields
if ~isfield(tendon1, 'e_start'), tendon1.e_start = []; end
if ~isfield(tendon1, 'e_mid'),   tendon1.e_mid = [];   end
if ~isfield(tendon1, 'e_end'),   tendon1.e_end = [];   end
if ~isfield(tendon1, 'drape_point'), tendon1.drape_point = []; end
if ~isfield(tendon1, 'cover_bot'),   tendon1.cover_bot = 2.0; end
if ~isfield(tendon1, 'cover_top'),   tendon1.cover_top = 2.0; end

% Bonding
tendon1.bonding.type = 'full';
tendon1.bonding.bonded_zones = [0, 1.0];

% Combine tendons
prestress.tendons = {tendon1};
prestress.losses = 0.0;

%% SECTION 5: NON-PRESTRESSED REINFORCEMENT
rebar_top.y_position = [];
rebar_top.x_positions = [];
rebar_top.areas = [];
rebar_top.type = 'compression';

rebar_bot.y_position = [];
rebar_bot.x_positions = [];
rebar_bot.areas = [];
rebar_bot.type = 'tension';

stirrups.Av = [];
stirrups.spacing = [];
stirrups.fy = [];

reinforcement.longitudinal = {rebar_top, rebar_bot};
reinforcement.stirrups = stirrups;

%% SECTION 6: APPLIED LOADS
loads.self_weight = [];
loads.concrete_density = 150/1728;

loads.distributed = [];
loads.point = [];

loads.support_type = 'simple';
loads.supports = [0, 100];

%% Calculate section properties if not provided
section = calculateSectionProperties(section);

%% Process tendon profiles
prestress = processTendonProfiles(prestress, beam, section);

%% Check tendon is within section boundaries
checkTendonBoundary(prestress, beam, section);

%% Display summary
displayInputSummary(beam, section, materials, prestress, reinforcement, loads);

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function section = calculateSectionProperties(section)
% Calculate section properties from polygon vertices using shoelace formula

v = section.vertices;
n = size(v, 1);
v = [v; v(1,:)];

A = 0;
for i = 1:n
    A = A + (v(i,1) * v(i+1,2) - v(i+1,1) * v(i,2));
end
A = abs(A) / 2;

Cx = 0; Cy = 0;
for i = 1:n
    factor = (v(i,1) * v(i+1,2) - v(i+1,1) * v(i,2));
    Cx = Cx + (v(i,1) + v(i+1,1)) * factor;
    Cy = Cy + (v(i,2) + v(i+1,2)) * factor;
end
Cx = Cx / (6 * A);
Cy = Cy / (6 * A);

Ix = 0; Iy = 0; Ixy = 0;
for i = 1:n
    xi  = v(i,1) - Cx;   yi  = v(i,2) - Cy;
    xi1 = v(i+1,1) - Cx; yi1 = v(i+1,2) - Cy;
    factor = (xi * yi1 - xi1 * yi);
    Ix  = Ix  + (yi^2 + yi*yi1 + yi1^2) * factor;
    Iy  = Iy  + (xi^2 + xi*xi1 + xi1^2) * factor;
    Ixy = Ixy + (xi*yi1 + 2*xi*yi + 2*xi1*yi1 + xi1*yi) * factor;
end
Ix = abs(Ix) / 12;
Iy = abs(Iy) / 12;

if isempty(section.A),  section.A  = A;  end
if isempty(section.Ix), section.Ix = Ix; end
if isempty(section.Iy), section.Iy = Iy; end
if isempty(section.xc), section.xc = Cx; end
if isempty(section.yc), section.yc = Cy; end

y_max = max(section.vertices(:,2));
y_min = min(section.vertices(:,2));
if isempty(section.yt), section.yt = y_max - section.yc; end
if isempty(section.yb), section.yb = section.yc - y_min; end

end

%% ========================================================================
function prestress = processTendonProfiles(prestress, beam, section)
% Generate eccentricity and y-coordinate profiles for each tendon
%
% Supported profile types:
%   straight, linear, parabolic, parabolic_y, parabolic_3pt,
%   parabolic_sym, harped, custom

x  = beam.x;
L  = beam.L;
yc = section.yc;

for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};
    
    switch tendon.profile_type
        
        case 'straight'
            % ---- Constant y-position ----
            y_tendon = ones(size(x)) * tendon.y_position;
            e = yc - y_tendon;
            tendon.e_start = e(1);
            tendon.e_mid   = e(round(end/2));
            tendon.e_end   = e(end);
            
        case 'parabolic_sym'
            % ---- Symmetric parabola from midspan ----
            %   y(x) = a * (x - L/2)^2 + y_mid
            %
            % The profile is symmetric about x = L/2.
            % At the ends (x=0 and x=L), the tendon rises to y_end.
            % y_end is determined by concrete cover unless overridden.
            
            y_mid_val = tendon.y_mid;
            
            % Determine end y-position from cover or override
            if isfield(tendon, 'y_end_override') && ~isempty(tendon.y_end_override)
                y_end_val = tendon.y_end_override;
            else
                y_end_val = tendon.cover_bot;  % y at ends = cover from bottom
            end
            
            % Calculate parabola coefficient
            %   y(0) = a*(L/2)^2 + y_mid = y_end
            %   => a = (y_end - y_mid) / (L/2)^2
            a_coeff = (y_end_val - y_mid_val) / (L/2)^2;
            
            y_tendon = a_coeff * (x - L/2).^2 + y_mid_val;
            e = yc - y_tendon;
            
            tendon.e_start = e(1);
            tendon.e_mid   = e(round(end/2));
            tendon.e_end   = e(end);
            tendon.y_coeffs_sym = [a_coeff, y_mid_val];  % store for reference
            
            fprintf('  Tendon %d (parabolic_sym):\n', i);
            fprintf('    y(x) = %.6f * (x - %.1f)^2 + %.2f\n', a_coeff, L/2, y_mid_val);
            fprintf('    y(0) = y(L) = %.2f in,  y(L/2) = %.2f in\n', y_tendon(1), y_mid_val);
            
        case 'parabolic_y'
            % ---- y(x) = a*x^2 + b*x + c (direct coefficients) ----
            a_c = tendon.y_coeffs(1);
            b_c = tendon.y_coeffs(2);
            c_c = tendon.y_coeffs(3);
            
            y_tendon = a_c * x.^2 + b_c * x + c_c;
            e = yc - y_tendon;
            
            tendon.e_start = e(1);
            tendon.e_mid   = e(round(end/2));
            tendon.e_end   = e(end);
            
            fprintf('  Tendon %d (parabolic_y): y = %.6f*x^2 + %.4f*x + %.2f\n', ...
                i, a_c, b_c, c_c);
            
        case 'parabolic_3pt'
            % ---- Fit parabola through y_start, y_mid, y_end ----
            y0 = tendon.y_start;
            ym = tendon.y_mid;
            yL = tendon.y_end;
            
            c_c = y0;
            a_c = 2 * (y0 - 2*ym + yL) / L^2;
            b_c = (yL - y0) / L - a_c * L;
            
            y_tendon = a_c * x.^2 + b_c * x + c_c;
            e = yc - y_tendon;
            
            tendon.e_start = e(1);
            tendon.e_mid   = e(round(end/2));
            tendon.e_end   = e(end);
            tendon.y_coeffs = [a_c, b_c, c_c];
            
            fprintf('  Tendon %d (parabolic_3pt): y = %.6f*x^2 + %.4f*x + %.2f\n', ...
                i, a_c, b_c, c_c);
            
        case 'linear'
            % ---- Linear eccentricity ----
            if isfield(tendon, 'y_position') && ~isempty(tendon.y_position)
                e_constant = yc - tendon.y_position;
                e = ones(size(x)) * e_constant;
                y_tendon = ones(size(x)) * tendon.y_position;
                tendon.e_start = e_constant;
                tendon.e_mid   = e_constant;
                tendon.e_end   = e_constant;
            else
                e_s = tendon.e_start;
                e_e = tendon.e_end;
                e = e_s + (e_e - e_s) * x / L;
                y_tendon = yc - e;
                tendon.e_mid = (e_s + e_e) / 2;
            end
            
        case 'parabolic'
            % ---- Parabolic eccentricity ----
            e_s = tendon.e_start;
            e_m = tendon.e_mid;
            e_e = tendon.e_end;
            
            a = 2 * (e_s - 2*e_m + e_e) / L^2;
            b = (e_e - e_s) / L - a * L;
            c = e_s;
            e = a * x.^2 + b * x + c;
            y_tendon = yc - e;
            
        case 'harped'
            % ---- Two linear segments at drape point ----
            e_s = tendon.e_start;
            e_m = tendon.e_mid;
            e_e = tendon.e_end;
            x_drape = tendon.drape_point * L;
            
            e = zeros(size(x));
            mask1 = x <= x_drape;
            mask2 = x > x_drape;
            e(mask1) = e_s + (e_m - e_s) * x(mask1) / x_drape;
            e(mask2) = e_m + (e_e - e_m) * (x(mask2) - x_drape) / (L - x_drape);
            y_tendon = yc - e;
            
        case 'custom'
            % ---- User-defined profile ----
            if isfield(tendon, 'e_profile')
                e = interp1(tendon.x_profile, tendon.e_profile, x, 'linear');
                y_tendon = yc - e;
                tendon.e_start = e(1);
                tendon.e_mid   = e(round(end/2));
                tendon.e_end   = e(end);
            else
                error('Custom profile requires e_profile and x_profile fields');
            end
            
        otherwise
            error('Unknown profile type: %s', tendon.profile_type);
    end
    
    % Store results back into tendon struct
    prestress.tendons{i} = tendon;
    prestress.tendons{i}.e = e;
    prestress.tendons{i}.y = y_tendon;
    prestress.tendons{i}.bonded = generateBondingMask(tendon.bonding, x, L);
end

end

%% ========================================================================
function checkTendonBoundary(prestress, beam, section)
% CHECK THAT ALL TENDONS REMAIN WITHIN THE CONCRETE SECTION
%
% For each tendon at each x-position, verify:
%   (y_bot + cover_bot) <= y_tendon <= (y_top - cover_top)
%
% Also checks that tendon y-coordinate is within the section height
% (i.e., not outside the polygon entirely).

x = beam.x;
y_bot = min(section.vertices(:,2));  % Bottom of section
y_top = max(section.vertices(:,2));  % Top of section

fprintf('\n--- Tendon Boundary Check ---\n');

for i = 1:length(prestress.tendons)
    t = prestress.tendons{i};
    y_tendon = t.y;
    
    % Get cover limits
    if isfield(t, 'cover_bot') && ~isempty(t.cover_bot)
        cov_bot = t.cover_bot;
    else
        cov_bot = 0;
    end
    if isfield(t, 'cover_top') && ~isempty(t.cover_top)
        cov_top = t.cover_top;
    else
        cov_top = 0;
    end
    
    y_min_allow = y_bot + cov_bot;
    y_max_allow = y_top - cov_top;
    
    has_violation = false;
    
    % --- Check bottom cover violation ---
    below_mask = y_tendon < y_min_allow;
    if any(below_mask)
        has_violation = true;
        idx_worst = find(y_tendon == min(y_tendon(below_mask)), 1);
        warning('TENDON %d: Violates BOTTOM cover at %d locations!', i, sum(below_mask));
        fprintf('  WARNING: Tendon %d below min y at %d points.\n', i, sum(below_mask));
        fprintf('    Worst: y = %.2f in at x = %.1f in (min allowed = %.2f in)\n', ...
            y_tendon(idx_worst), x(idx_worst), y_min_allow);
        fprintf('    Required cover = %.2f in, actual = %.2f in\n', ...
            cov_bot, y_tendon(idx_worst) - y_bot);
    end
    
    % --- Check top cover violation ---
    above_mask = y_tendon > y_max_allow;
    if any(above_mask)
        has_violation = true;
        idx_worst = find(y_tendon == max(y_tendon(above_mask)), 1);
        warning('TENDON %d: Violates TOP cover at %d locations!', i, sum(above_mask));
        fprintf('  WARNING: Tendon %d above max y at %d points.\n', i, sum(above_mask));
        fprintf('    Worst: y = %.2f in at x = %.1f in (max allowed = %.2f in)\n', ...
            y_tendon(idx_worst), x(idx_worst), y_max_allow);
        fprintf('    Required cover = %.2f in, actual = %.2f in\n', ...
            cov_top, y_top - y_tendon(idx_worst));
    end
    
    % --- Check tendon is inside section height entirely ---
    outside_bot = y_tendon < y_bot;
    outside_top = y_tendon > y_top;
    if any(outside_bot) || any(outside_top)
        has_violation = true;
        n_out = sum(outside_bot) + sum(outside_top);
        warning('TENDON %d: OUTSIDE the section at %d locations!', i, n_out);
        fprintf('  WARNING: Tendon %d is completely outside the concrete section!\n', i);
        if any(outside_bot)
            idx = find(outside_bot);
            fprintf('    Below section at x = [');
            fprintf('%.1f ', x(idx(1:min(5,end))));
            if length(idx) > 5, fprintf('...'); end
            fprintf('] in\n');
        end
        if any(outside_top)
            idx = find(outside_top);
            fprintf('    Above section at x = [');
            fprintf('%.1f ', x(idx(1:min(5,end))));
            if length(idx) > 5, fprintf('...'); end
            fprintf('] in\n');
        end
    end
    
    % --- Check tendon is within section width at its y-level ---
    % For I-sections, the web is narrower than flanges
    outside_width = false(size(x));
    for j = 1:length(x)
        y_t = y_tendon(j);
        if y_t >= y_bot && y_t <= y_top
            half_w = getSectionHalfWidth(section.vertices, y_t);
            if isnan(half_w) || half_w < 0.5  % tendon needs some minimum clearance
                outside_width(j) = true;
            end
        end
    end
    if any(outside_width)
        has_violation = true;
        idx_out = find(outside_width);
        warning('TENDON %d: At y-levels where section is too narrow at %d locations!', i, length(idx_out));
        fprintf('  WARNING: Tendon %d at y-levels where section width is insufficient.\n', i);
        fprintf('    Locations: x = [');
        fprintf('%.1f ', x(idx_out(1:min(5,end))));
        if length(idx_out) > 5, fprintf('...'); end
        fprintf('] in\n');
    end
    
    % --- All OK ---
    if ~has_violation
        fprintf('  Tendon %d: OK - within section with adequate cover\n', i);
        fprintf('    y range: [%.2f, %.2f] in,  allowed: [%.2f, %.2f] in\n', ...
            min(y_tendon), max(y_tendon), y_min_allow, y_max_allow);
    end
end

fprintf('--- End Boundary Check ---\n\n');

end

%% ========================================================================
function half_w = getSectionHalfWidth(vertices, y_query)
% Find the half-width of the section at a given y-coordinate
% Returns NaN if y_query is outside the section

y_min = min(vertices(:,2));
y_max = max(vertices(:,2));

if y_query < y_min || y_query > y_max
    half_w = NaN;
    return;
end

% Find all x-coordinates where horizontal line y=y_query intersects polygon
n = size(vertices, 1);
x_intersections = [];

for i = 1:n
    j = mod(i, n) + 1;
    y1 = vertices(i, 2);
    y2 = vertices(j, 2);
    
    if (y1 <= y_query && y2 > y_query) || (y2 <= y_query && y1 > y_query)
        t = (y_query - y1) / (y2 - y1);
        x_int = vertices(i, 1) + t * (vertices(j, 1) - vertices(i, 1));
        x_intersections = [x_intersections, x_int]; %#ok<AGROW>
    end
end

if isempty(x_intersections)
    half_w = NaN;
else
    half_w = (max(x_intersections) - min(x_intersections)) / 2;
end

end

%% ========================================================================
function bonded = generateBondingMask(bonding, x, L)

bonded = false(size(x));

switch bonding.type
    case 'full'
        bonded(:) = true;
    case 'unbonded'
        bonded(:) = false;
    case 'partial'
        for i = 1:size(bonding.bonded_zones, 1)
            x_start = bonding.bonded_zones(i, 1) * L;
            x_end   = bonding.bonded_zones(i, 2) * L;
            bonded = bonded | (x >= x_start & x <= x_end);
        end
end

end

%% ========================================================================
function displayInputSummary(beam, section, materials, prestress, reinforcement, loads)

fprintf('\n========================================\n');
fprintf('   HOMEWORK 1 - PRESTRESSED BEAM INPUT\n');
fprintf('========================================\n');

fprintf('\nBEAM GEOMETRY:\n');
fprintf('  Span length: %.1f in (%.2f ft)\n', beam.L, beam.L/12);
fprintf('  Analysis segments: %d\n', beam.num_segments);

fprintf('\nSECTION PROPERTIES:\n');
fprintf('  Area (A): %.2f in^2\n', section.A);
fprintf('  Moment of inertia (Ix): %.2f in^4\n', section.Ix);
fprintf('  Centroid from bottom (yc): %.2f in\n', section.yc);
fprintf('  Centroid from left (xc): %.2f in\n', section.xc);
fprintf('  Distance to top fiber (yt): %.2f in\n', section.yt);
fprintf('  Distance to bottom fiber (yb): %.2f in\n', section.yb);

fprintf('\nMATERIAL PROPERTIES:\n');
fprintf('  Concrete: f''c = %.1f ksi (%.0f psi), Ec = %.0f ksi\n', ...
    materials.fc, materials.fc*1000, materials.Ec);
fprintf('  Prestressing steel: fpu = %.0f ksi, Eps = %.0f ksi\n', ...
    materials.fpu, materials.Eps);
fprintf('  Mild steel: fy = %.0f ksi, Es = %.0f ksi\n', ...
    materials.fy, materials.Es);

fprintf('\nPRESTRESSING TENDONS:\n');
total_Aps = 0;
total_P = 0;
for i = 1:length(prestress.tendons)
    t = prestress.tendons{i};
    P_initial = t.Aps * t.fpi;
    fprintf('  Tendon %d:\n', i);
    fprintf('    Aps = %.3f in^2\n', t.Aps);
    fprintf('    fpi = %.1f ksi\n', t.fpi);
    fprintf('    Initial force (P) = %.1f kips\n', P_initial);
    fprintf('    Profile type: %s\n', t.profile_type);
    
    switch t.profile_type
        case 'straight'
            fprintf('    y = %.2f in (constant)\n', t.y_position);
            fprintf('    e = %.2f in (constant)\n', t.e_start);
        case 'parabolic_sym'
            fprintf('    y(x) = %.6f*(x - %.1f)^2 + %.2f   [symmetric]\n', ...
                t.y_coeffs_sym(1), beam.L/2, t.y_coeffs_sym(2));
            fprintf('    y(0) = y(L) = %.2f in,  y(L/2) = %.2f in\n', ...
                t.y(1), t.y(round(end/2)));
            fprintf('    e: %.2f (end) to %.2f (mid) in\n', t.e_start, t.e_mid);
        case {'parabolic_y', 'parabolic_3pt'}
            if isfield(t, 'y_coeffs')
                fprintf('    y(x) = %.6f*x^2 + %.4f*x + %.2f\n', ...
                    t.y_coeffs(1), t.y_coeffs(2), t.y_coeffs(3));
            end
            fprintf('    y: %.2f (start) -> %.2f (mid) -> %.2f (end) in\n', ...
                t.y(1), t.y(round(end/2)), t.y(end));
            fprintf('    e: %.2f (start) -> %.2f (mid) -> %.2f (end) in\n', ...
                t.e_start, t.e_mid, t.e_end);
        otherwise
            fprintf('    e: %.2f (start) -> %.2f (mid) -> %.2f (end) in\n', ...
                t.e_start, t.e_mid, t.e_end);
    end
    
    % Cover info
    if isfield(t, 'cover_bot')
        fprintf('    Cover: bot = %.1f in, top = %.1f in\n', t.cover_bot, t.cover_top);
    end
    fprintf('    Bonding: %s\n', t.bonding.type);
    
    total_Aps = total_Aps + t.Aps;
    total_P = total_P + P_initial;
end
fprintf('  Total Aps: %.3f in^2\n', total_Aps);
fprintf('  Total P (initial): %.1f kips\n', total_P);
fprintf('  Prestress losses: %.1f%%\n', prestress.losses * 100);
fprintf('  Total P (after losses): %.1f kips\n', total_P * (1 - prestress.losses));

fprintf('\nNON-PRESTRESSED REINFORCEMENT:\n');
has_rebar = false;
for i = 1:length(reinforcement.longitudinal)
    r = reinforcement.longitudinal{i};
    if ~isempty(r.y_position) && ~isempty(r.areas)
        fprintf('  %s: y = %.1f in, As = %.2f in^2\n', ...
            r.type, r.y_position, sum(r.areas));
        has_rebar = true;
    end
end
if ~has_rebar
    fprintf('  None specified\n');
end

fprintf('\nAPPLIED LOADS:\n');
if ~isempty(loads.distributed)
    fprintf('  Distributed loads: %d segments\n', size(loads.distributed, 1));
else
    fprintf('  Distributed loads: None\n');
end
if ~isempty(loads.point)
    fprintf('  Point loads: %d\n', size(loads.point, 1));
else
    fprintf('  Point loads: None\n');
end
fprintf('  Support type: %s\n', loads.support_type);

fprintf('\n========================================\n');

end
