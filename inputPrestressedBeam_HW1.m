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
Aps_required = 400 / fpi;       % Area required for 400 kips (in^2)

% Tendon 1 - Straight profile at y=6 (fully bonded)
tendon1.Aps = Aps_required;     % Area of prestressing steel (in^2)
tendon1.fpi = fpi;              % Initial prestress (ksi)
tendon1.profile_type = 'linear'; % Straight tendon
tendon1.x_position = 0;         % x-coordinate of tendon (in)
tendon1.y_position = 6;         % y-coordinate of tendon from bottom (in)

% Eccentricity will be calculated as: e = yc - y_tendon
% For straight profile, e is constant along the beam
tendon1.e_start = [];           % Will be calculated from y_position
tendon1.e_mid = [];             % Will be calculated from y_position
tendon1.e_end = [];             % Will be calculated from y_position
tendon1.drape_point = [];       % Not applicable for straight profile

% Bonding information for Tendon 1 - FULLY BONDED
tendon1.bonding.type = 'full';
tendon1.bonding.bonded_zones = [0, 1.0];  % [start, end] as fraction of L

% Combine tendons
prestress.tendons = {tendon1};
prestress.losses = 0.15;        % Total prestress losses as fraction (15%)

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

x = beam.x;
L = beam.L;
yc = section.yc;

for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};
    
    % For straight tendon at fixed y-position
    if isfield(tendon, 'y_position') && ~isempty(tendon.y_position)
        % Calculate constant eccentricity: e = yc - y_tendon
        e_constant = yc - tendon.y_position;
        
        % Set all eccentricities to constant value
        tendon.e_start = e_constant;
        tendon.e_mid = e_constant;
        tendon.e_end = e_constant;
        
        % Create constant eccentricity profile
        e = ones(size(x)) * e_constant;
        
        % Store y-coordinate
        y_tendon = ones(size(x)) * tendon.y_position;
    else
        % Use provided eccentricity values
        e_s = tendon.e_start;
        e_m = tendon.e_mid;
        e_e = tendon.e_end;
        
        switch tendon.profile_type
            case 'linear'
                % Linear variation from start to end
                e = e_s + (e_e - e_s) * x / L;
                
            case 'parabolic'
                % Parabolic profile: e(x) = ax^2 + bx + c
                % Boundary conditions: e(0)=e_s, e(L/2)=e_m, e(L)=e_e
                a = 2 * (e_s - 2*e_m + e_e) / L^2;
                b = (e_e - e_s) / L - a * L;
                c = e_s;
                e = a * x.^2 + b * x + c;
                
            case 'harped'
                % Harped profile: linear segments meeting at drape point
                x_drape = tendon.drape_point * L;
                e = zeros(size(x));
                mask1 = x <= x_drape;
                mask2 = x > x_drape;
                e(mask1) = e_s + (e_m - e_s) * x(mask1) / x_drape;
                e(mask2) = e_m + (e_e - e_m) * (x(mask2) - x_drape) / (L - x_drape);
                
            case 'custom'
                % User-defined profile (requires tendon.e_profile field)
                if isfield(tendon, 'e_profile')
                    e = interp1(tendon.x_profile, tendon.e_profile, x, 'linear');
                else
                    error('Custom profile requires e_profile and x_profile fields');
                end
                
            otherwise
                error('Unknown profile type: %s', tendon.profile_type);
        end
        
        % Calculate tendon y-coordinate from bottom
        y_tendon = yc - e;
    end
    
    % Store eccentricity profile (positive = below centroid)
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
    if isfield(t, 'y_position')
        fprintf('    Location: (x=%.1f, y=%.1f) in\n', t.x_position, t.y_position);
        fprintf('    Eccentricity: e = %.2f in (constant)\n', t.e_start);
    else
        fprintf('    Profile: %s\n', t.profile_type);
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
