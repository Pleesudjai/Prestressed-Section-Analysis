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
%
% Outputs:
%   beam - Structure containing beam geometry and span information
%   section - Structure containing cross-section geometry (polygonal)
%   materials - Structure containing material properties
%   prestress - Structure containing prestressing tendon information
%   reinforcement - Structure containing non-prestressed reinforcement
%   loads - Structure containing applied loads

%% SECTION 1: BEAM GEOMETRY
beam.L = 100;                   % Total beam length (in)
beam.num_segments = 100;        % Number of segments for analysis
beam.x = linspace(0, beam.L, beam.num_segments + 1);  % x-coordinates along beam

%% SECTION 2: CROSS-SECTION GEOMETRY
% Define vertices counterclockwise [x, y] - I-section centered at origin
% Coordinates from homework 1 specification
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

% Section properties (will be calculated automatically if empty)
section.A = [];                 % Cross-sectional area (in^2)
section.Ix = [];                % Moment of inertia about x-axis (in^4)
section.Iy = [];                % Moment of inertia about y-axis (in^4)
section.yc = [];                % Centroid y-coordinate from bottom (in)
section.xc = [];                % Centroid x-coordinate from left (in)
section.yt = [];                % Distance from centroid to top fiber (in)
section.yb = [];                % Distance from centroid to bottom fiber (in)

%% SECTION 3: MATERIAL PROPERTIES
% Concrete properties - from homework specification
materials.fc = 6.0;             % Concrete compressive strength (ksi) = 6000 psi
materials.Ec = 4700;            % Concrete modulus (ksi) - given directly
materials.fci = 4.5;            % Concrete strength at transfer (ksi)
materials.fr = 7.5 * sqrt(materials.fc * 1000) / 1000;  % Modulus of rupture (ksi)

% Prestressing steel properties
materials.fpu = 270;            % Ultimate strength of prestressing steel (ksi)
materials.fpy = 0.9 * materials.fpu;  % Yield strength (ksi)
materials.Eps = 28500;          % Modulus of elasticity of prestressing steel (ksi)

% Non-prestressed steel properties
materials.fy = 60;              % Yield strength of mild steel (ksi)
materials.Es = 29000;           % Modulus of elasticity of mild steel (ksi)

%% SECTION 4: PRESTRESSING TENDONS
% Single tendon at (x=0, y=6) with P=400 kips, fully bonded
% Need to calculate area: Aps = P / fpi

% Assume initial prestress at 75% of fpu
fpi = 0.75 * materials.fpu;     % Initial prestress (ksi)
Aps_required = 400 / fpi ;       % Area required for 400 kips (in^2)

% =====================================================================
% TENDON PROFILE TYPES - Choose one of the following:
%
%   'straight'       - Constant y-position along beam
%                      Requires: y_position
%
%   'linear'         - Linear y variation from start to end
%                      Requires: e_start, e_end (eccentricities)
%
%   'parabolic'      - Parabolic eccentricity: e(x) = ax^2 + bx + c
%                      Requires: e_start, e_mid, e_end
%
%   'parabolic_y'    - Parabolic y-coordinate: y(x) = a*x^2 + b*x + c
%                      Requires: y_coeffs = [a, b, c]
%                      Example: y = -0.001*x^2 + 0.1*x + 6
%
%   'parabolic_3pt'  - Parabolic y from 3 control points
%                      Requires: y_start, y_mid, y_end
%                      Fits y(x) = a*x^2 + b*x + c through the 3 points
%
%   'harped'         - Two linear segments meeting at drape point
%                      Requires: e_start, e_mid, e_end, drape_point
%
%   'custom'         - User-defined profile
%                      Requires: x_profile, e_profile (vectors)
% =====================================================================

% Tendon 1 - Choose profile type here
tendon1.Aps = Aps_required;     % Area of prestressing steel (in^2)
tendon1.fpi = fpi;              % Initial prestress (ksi)

% --- ACTIVE PROFILE (uncomment ONE block below) ---

% % OPTION A: Straight tendon at constant y
% tendon1.profile_type = 'straight';
% tendon1.y_position = 6;         % y-coordinate of tendon from bottom (in)

% OPTION B: Parabolic y-coordinate using coefficients y = a*x^2 + b*x + c
tendon1.profile_type = 'parabolic_y';
tendon1.y_coeffs = [1, 0, 0];  % [a, b, c] -> y = a*x^2 + b*x + c
  % This example: y(0)=6, y(50)=10, y(100)=6 in

% % OPTION C: Parabolic y from 3 control points (start, midspan, end)
% tendon1.profile_type = 'parabolic_3pt';
% tendon1.y_start = 6;           % y at x = 0
% tendon1.y_mid   = 10;          % y at x = L/2
% tendon1.y_end   = 6;           % y at x = L

% % OPTION D: Parabolic eccentricity e(x) = ax^2 + bx + c
% tendon1.profile_type = 'parabolic';
% tendon1.e_start = 30;          % Eccentricity at x = 0
% tendon1.e_mid   = 35;          % Eccentricity at x = L/2
% tendon1.e_end   = 30;          % Eccentricity at x = L

% % OPTION E: Harped profile
% tendon1.profile_type = 'harped';
% tendon1.e_start = 10;
% tendon1.e_mid   = 35;
% tendon1.e_end   = 10;
% tendon1.drape_point = 0.5;     % Fraction of L where drape occurs

% Common fields (defaults for unused options)
if ~isfield(tendon1, 'e_start'), tendon1.e_start = []; end
if ~isfield(tendon1, 'e_mid'),   tendon1.e_mid = [];   end
if ~isfield(tendon1, 'e_end'),   tendon1.e_end = [];   end
if ~isfield(tendon1, 'drape_point'), tendon1.drape_point = []; end

% Bonding information for Tendon 1 - FULLY BONDED
tendon1.bonding.type = 'full';
tendon1.bonding.bonded_zones = [0, 1.0];  % [start, end] as fraction of L

% Combine tendons
prestress.tendons = {tendon1};
prestress.losses = 0.0; %0.15;        % Total prestress losses as fraction (15%)

%% SECTION 5: NON-PRESTRESSED REINFORCEMENT
% Define mild steel reinforcement (if any)
% For now, assume no additional non-prestressed reinforcement

% Top reinforcement (compression steel) - empty for now
rebar_top.y_position = [];
rebar_top.x_positions = [];
rebar_top.areas = [];
rebar_top.type = 'compression';

% Bottom reinforcement (tension steel) - empty for now
rebar_bot.y_position = [];
rebar_bot.x_positions = [];
rebar_bot.areas = [];
rebar_bot.type = 'tension';

% Shear reinforcement (stirrups) - empty for now
stirrups.Av = [];
stirrups.spacing = [];
stirrups.fy = [];

% Combine reinforcement
reinforcement.longitudinal = {rebar_top, rebar_bot};
reinforcement.stirrups = stirrups;

%% SECTION 6: APPLIED LOADS
% Self-weight (calculated automatically if empty)
loads.self_weight = [];         % Self-weight (kip/in), leave empty for auto-calc
loads.concrete_density = 150/1728;  % Concrete density (kip/in^3)

% Distributed loads [start_x, end_x, w_start, w_end] (kip/in)
% Each row defines a load segment: [x_start, x_end, w_at_start, w_at_end]
% For homework, may need to add specific loads later
loads.distributed = [
    % 0,      600,    -0.05,  -0.05;   % Example: Uniform dead load (negative = downward)
];

% Point loads [x_position, P, M] (kips and kip-in)
% P: positive = upward, M: positive = counterclockwise
loads.point = [
    % 150,    -20,    0;      % Example: Point load
];

% Support conditions: 'simple', 'cantilever', 'continuous', 'custom'
loads.support_type = 'simple';

% For simple supports: [x_left, x_right]
loads.supports = [0, 100];

%% Calculate section properties if not provided
section = calculateSectionProperties(section);

%% Process tendon profiles - update for straight tendon at y=6
prestress = processTendonProfiles(prestress, beam, section);

%% Display summary
displayInputSummary(beam, section, materials, prestress, reinforcement, loads);

end

%% HELPER FUNCTIONS

function section = calculateSectionProperties(section)
% Calculate section properties from polygon vertices using shoelace formula

v = section.vertices;
n = size(v, 1);

% Close the polygon
v = [v; v(1,:)];

% Calculate area using shoelace formula
A = 0;
for i = 1:n
    A = A + (v(i,1) * v(i+1,2) - v(i+1,1) * v(i,2));
end
A = abs(A) / 2;

% Calculate centroid
Cx = 0; Cy = 0;
for i = 1:n
    factor = (v(i,1) * v(i+1,2) - v(i+1,1) * v(i,2));
    Cx = Cx + (v(i,1) + v(i+1,1)) * factor;
    Cy = Cy + (v(i,2) + v(i+1,2)) * factor;
end
Cx = Cx / (6 * A);
Cy = Cy / (6 * A);

% Calculate moments of inertia about centroidal axes
Ix = 0; Iy = 0; Ixy = 0;
for i = 1:n
    xi = v(i,1) - Cx;
    yi = v(i,2) - Cy;
    xi1 = v(i+1,1) - Cx;
    yi1 = v(i+1,2) - Cy;
    
    factor = (xi * yi1 - xi1 * yi);
    Ix = Ix + (yi^2 + yi*yi1 + yi1^2) * factor;
    Iy = Iy + (xi^2 + xi*xi1 + xi1^2) * factor;
    Ixy = Ixy + (xi*yi1 + 2*xi*yi + 2*xi1*yi1 + xi1*yi) * factor;
end
Ix = abs(Ix) / 12;
Iy = abs(Iy) / 12;

% Assign calculated values if not provided
if isempty(section.A), section.A = A; end
if isempty(section.Ix), section.Ix = Ix; end
if isempty(section.Iy), section.Iy = Iy; end
if isempty(section.xc), section.xc = Cx; end
if isempty(section.yc), section.yc = Cy; end

% Calculate distances to extreme fibers
y_max = max(section.vertices(:,2));
y_min = min(section.vertices(:,2));
if isempty(section.yt), section.yt = y_max - section.yc; end
if isempty(section.yb), section.yb = section.yc - y_min; end

end

function prestress = processTendonProfiles(prestress, beam, section)
% Generate eccentricity profile along beam length for each tendon
%
% Supported profile_type values:
%   'straight'      - Constant y_position
%   'linear'        - Linear eccentricity from e_start to e_end
%   'parabolic'     - Parabolic eccentricity from e_start, e_mid, e_end
%   'parabolic_y'   - Parabolic y(x) = a*x^2 + b*x + c via y_coeffs
%   'parabolic_3pt' - Parabolic y fitted through y_start, y_mid, y_end
%   'harped'        - Two linear segments with drape point
%   'custom'        - User-defined e_profile and x_profile vectors

x = beam.x;
L = beam.L;
yc = section.yc;

for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};
    
    switch tendon.profile_type
        
        case 'straight'
            % Constant y-position along beam
            y_tendon = ones(size(x)) * tendon.y_position;
            e = yc - y_tendon;
            
            % Store eccentricity summary values
            tendon.e_start = e(1);
            tendon.e_mid = e(round(end/2));
            tendon.e_end = e(end);
            
        case 'parabolic_y'
            % y(x) = a*x^2 + b*x + c  (direct coefficient input)
            a_coeff = tendon.y_coeffs(1);
            b_coeff = tendon.y_coeffs(2);
            c_coeff = tendon.y_coeffs(3);
            
            y_tendon = a_coeff * x.^2 + b_coeff * x + c_coeff;
            e = yc - y_tendon;
            
            tendon.e_start = e(1);
            tendon.e_mid = e(round(end/2));
            tendon.e_end = e(end);
            
            fprintf('  Tendon %d (parabolic_y): y = %.6f*x^2 + %.4f*x + %.2f\n', ...
                i, a_coeff, b_coeff, c_coeff);
            fprintf('    y(0)=%.2f, y(L/2)=%.2f, y(L)=%.2f in\n', ...
                y_tendon(1), y_tendon(round(end/2)), y_tendon(end));
            
        case 'parabolic_3pt'
            % Fit y(x) = a*x^2 + b*x + c through 3 points:
            %   y(0) = y_start,  y(L/2) = y_mid,  y(L) = y_end
            y0 = tendon.y_start;
            ym = tendon.y_mid;
            yL = tendon.y_end;
            
            % Solve: c = y0
            %        a*(L/2)^2 + b*(L/2) + c = ym
            %        a*L^2     + b*L     + c = yL
            c_coeff = y0;
            % From the two remaining equations:
            %   a*L^2/4 + b*L/2 = ym - y0
            %   a*L^2   + b*L   = yL - y0
            % Multiply first by 2: a*L^2/2 + b*L = 2*(ym - y0)
            % Subtract from second: a*L^2/2 = (yL - y0) - 2*(ym - y0)
            a_coeff = 2 * (y0 - 2*ym + yL) / L^2;
            b_coeff = (yL - y0) / L - a_coeff * L;
            
            y_tendon = a_coeff * x.^2 + b_coeff * x + c_coeff;
            e = yc - y_tendon;
            
            tendon.e_start = e(1);
            tendon.e_mid = e(round(end/2));
            tendon.e_end = e(end);
            
            tendon.y_coeffs = [a_coeff, b_coeff, c_coeff];
            
            fprintf('  Tendon %d (parabolic_3pt): y = %.6f*x^2 + %.4f*x + %.2f\n', ...
                i, a_coeff, b_coeff, c_coeff);
            fprintf('    y(0)=%.2f, y(L/2)=%.2f, y(L)=%.2f in\n', ...
                y_tendon(1), y_tendon(round(end/2)), y_tendon(end));
            
        case 'linear'
            % Linear eccentricity variation
            if isfield(tendon, 'y_position') && ~isempty(tendon.y_position)
                % Straight tendon (backward compatibility)
                e_constant = yc - tendon.y_position;
                e = ones(size(x)) * e_constant;
                y_tendon = ones(size(x)) * tendon.y_position;
                tendon.e_start = e_constant;
                tendon.e_mid = e_constant;
                tendon.e_end = e_constant;
            else
                e_s = tendon.e_start;
                e_e = tendon.e_end;
                e = e_s + (e_e - e_s) * x / L;
                y_tendon = yc - e;
                tendon.e_mid = (e_s + e_e) / 2;
            end
            
        case 'parabolic'
            % Parabolic eccentricity: e(x) = ax^2 + bx + c
            % Boundary conditions: e(0)=e_s, e(L/2)=e_m, e(L)=e_e
            e_s = tendon.e_start;
            e_m = tendon.e_mid;
            e_e = tendon.e_end;
            
            a = 2 * (e_s - 2*e_m + e_e) / L^2;
            b = (e_e - e_s) / L - a * L;
            c = e_s;
            e = a * x.^2 + b * x + c;
            y_tendon = yc - e;
            
        case 'harped'
            % Harped profile: linear segments meeting at drape point
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
            % User-defined profile (requires tendon.e_profile and x_profile)
            if isfield(tendon, 'e_profile')
                e = interp1(tendon.x_profile, tendon.e_profile, x, 'linear');
                y_tendon = yc - e;
                tendon.e_start = e(1);
                tendon.e_mid = e(round(end/2));
                tendon.e_end = e(end);
            else
                error('Custom profile requires e_profile and x_profile fields');
            end
            
        otherwise
            error('Unknown profile type: %s', tendon.profile_type);
    end
    
    % Store eccentricity profile (positive = below centroid)
    prestress.tendons{i} = tendon;
    prestress.tendons{i}.e = e;
    
    % Store tendon y-coordinate
    prestress.tendons{i}.y = y_tendon;
    
    % Generate bonding mask
    prestress.tendons{i}.bonded = generateBondingMask(tendon.bonding, x, L);
end

end

function bonded = generateBondingMask(bonding, x, L)
% Generate boolean mask indicating bonded regions

bonded = false(size(x));

switch bonding.type
    case 'full'
        bonded(:) = true;
        
    case 'unbonded'
        bonded(:) = false;
        
    case 'partial'
        % bonded_zones is [n x 2] matrix of [start, end] fractions
        for i = 1:size(bonding.bonded_zones, 1)
            x_start = bonding.bonded_zones(i, 1) * L;
            x_end = bonding.bonded_zones(i, 2) * L;
            bonded = bonded | (x >= x_start & x <= x_end);
        end
end

end

function displayInputSummary(beam, section, materials, prestress, reinforcement, loads)
% Display summary of input data

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
            fprintf('    Location: y = %.1f in (constant)\n', t.y_position);
            fprintf('    Eccentricity: e = %.2f in (constant)\n', t.e_start);
        case {'parabolic_y', 'parabolic_3pt'}
            if isfield(t, 'y_coeffs')
                fprintf('    y(x) = %.6f*x^2 + %.4f*x + %.2f\n', ...
                    t.y_coeffs(1), t.y_coeffs(2), t.y_coeffs(3));
            end
            fprintf('    y: %.2f in (start) to %.2f in (mid) to %.2f in (end)\n', ...
                t.y(1), t.y(round(end/2)), t.y(end));
            fprintf('    Eccentricity: %.2f in (start) to %.2f in (mid) to %.2f in (end)\n', ...
                t.e_start, t.e_mid, t.e_end);
        otherwise
            fprintf('    Eccentricity: %.2f in (start) to %.2f in (mid) to %.2f in (end)\n', ...
                t.e_start, t.e_mid, t.e_end);
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
