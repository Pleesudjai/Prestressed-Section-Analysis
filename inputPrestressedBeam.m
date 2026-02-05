function [beam, section, materials, prestress, reinforcement, loads] = inputPrestressedBeam()
% INPUTPRESTRESSEDBEAM - Input data for prestressed concrete beam analysis
% Creates axial force, shear force, and bending moment diagrams
%
% Outputs:
%   beam - Structure containing beam geometry and span information
%   section - Structure containing cross-section geometry (polygonal)
%   materials - Structure containing material properties
%   prestress - Structure containing prestressing tendon information
%   reinforcement - Structure containing non-prestressed reinforcement
%   loads - Structure containing applied loads

%% SECTION 1: BEAM GEOMETRY
beam.L = 600;                   % Total beam length (in)
beam.num_segments = 100;        % Number of segments for analysis
beam.x = linspace(0, beam.L, beam.num_segments + 1);  % x-coordinates along beam

%% SECTION 2: CROSS-SECTION GEOMETRY
% Define vertices counterclockwise [x, y] - origin at bottom-left
% Example: I-section
section.vertices = [
    0,    0;        % Bottom-left of bottom flange
    24,   0;        % Bottom-right of bottom flange
    24,   6;        % Top-right of bottom flange
    15,   6;        % Bottom-right of web
    15,   30;       % Top-right of web
    24,   30;       % Bottom-right of top flange
    24,   36;       % Top-right of top flange
    0,    36;       % Top-left of top flange
    0,    30;       % Bottom-left of top flange
    9,    30;       % Top-left of web
    9,    6;        % Bottom-left of web
    0,    6;        % Top-left of bottom flange
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
% Concrete properties
materials.fc = 6;               % Concrete compressive strength (ksi)
materials.Ec = 57 * sqrt(materials.fc * 1000);  % Concrete modulus (ksi)
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
% Define tendon profiles and properties
% Each tendon has: area, initial stress, and profile (eccentricity vs x)

% Tendon 1 - Parabolic profile
tendon1.Aps = 2.0;              % Area of prestressing steel (in^2)
tendon1.fpi = 0.75 * materials.fpu;  % Initial prestress (ksi)
tendon1.profile_type = 'parabolic';  % 'linear', 'parabolic', 'harped', 'custom'
tendon1.e_start = 4;            % Eccentricity at start (in, positive below centroid)
tendon1.e_mid = 12;             % Eccentricity at midspan (in)
tendon1.e_end = 4;              % Eccentricity at end (in)
tendon1.drape_point = 0.5;      % Drape point as fraction of span (for harped)

% Bonding information for Tendon 1
tendon1.bonding.type = 'partial';  % 'full', 'unbonded', 'partial'
tendon1.bonding.bonded_zones = [0, 0.1; 0.4, 0.6; 0.9, 1.0];  % [start, end] as fraction of L

% Tendon 2 - Straight profile (optional, can be commented out)
tendon2.Aps = 1.0;              % Area of prestressing steel (in^2)
tendon2.fpi = 0.70 * materials.fpu;  % Initial prestress (ksi)
tendon2.profile_type = 'linear';
tendon2.e_start = 6;            % Eccentricity at start (in)
tendon2.e_mid = 6;              % Eccentricity at midspan (in)
tendon2.e_end = 6;              % Eccentricity at end (in)
tendon2.drape_point = [];

% Bonding information for Tendon 2
tendon2.bonding.type = 'full';     % Fully bonded
tendon2.bonding.bonded_zones = [0, 1.0];

% Combine tendons
prestress.tendons = {tendon1, tendon2};
prestress.losses = 0.0 %0.15;        % Total prestress losses as fraction (15%)

%% SECTION 5: NON-PRESTRESSED REINFORCEMENT
% Define mild steel reinforcement

% Top reinforcement (compression steel)
rebar_top.y_position = 34;      % y-coordinate from bottom (in)
rebar_top.x_positions = [3, 12, 21];  % x-coordinates of bars (in)
rebar_top.areas = [0.44, 0.44, 0.44];  % Bar areas (in^2) - #6 bars
rebar_top.type = 'compression';

% Bottom reinforcement (tension steel)
rebar_bot.y_position = 2;       % y-coordinate from bottom (in)
rebar_bot.x_positions = [3, 8, 16, 21];  % x-coordinates of bars (in)
rebar_bot.areas = [0.60, 0.60, 0.60, 0.60];  % Bar areas (in^2) - #7 bars
rebar_bot.type = 'tension';

% Shear reinforcement (stirrups)
stirrups.Av = 0.22;             % Area of two legs (in^2) - #3 stirrups
stirrups.spacing = 6;           % Spacing (in)
stirrups.fy = 60;               % Yield strength (ksi)

% Combine reinforcement
reinforcement.longitudinal = {rebar_top, rebar_bot};
reinforcement.stirrups = stirrups;

%% SECTION 6: APPLIED LOADS
% Self-weight (calculated automatically if empty)
loads.self_weight = [];         % Self-weight (kip/in), leave empty for auto-calc
loads.concrete_density = 150/1728;  % Concrete density (kip/in^3)

% Distributed loads [start_x, end_x, w_start, w_end] (kip/in)
% Each row defines a load segment: [x_start, x_end, w_at_start, w_at_end]
loads.distributed = [
    0,      600,    -0.05,  -0.05;   % Uniform dead load (negative = downward)
    100,    500,    -0.03,  -0.03;   % Partial uniform live load
];

% Point loads [x_position, P, M] (kips and kip-in)
% P: positive = upward, M: positive = counterclockwise
loads.point = [
    150,    -20,    0;      % Point load at x=150
    450,    -20,    0;      % Point load at x=450
];

% Support conditions: 'simple', 'cantilever', 'continuous', 'custom'
loads.support_type = 'simple';

% For simple supports: [x_left, x_right]
loads.supports = [0, 600];

% For custom supports: define reaction locations and types
% loads.custom_supports = [
%     0,      'pin';      % Pin at left
%     300,    'roller';   % Roller at midspan
%     600,    'pin';      % Pin at right
% ];

%% Calculate section properties if not provided
section = calculateSectionProperties(section);

%% Process tendon profiles
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
    
    % Store eccentricity profile (positive = below centroid)
    prestress.tendons{i}.e = e;
    
    % Calculate tendon y-coordinate from bottom
    prestress.tendons{i}.y = yc - e;
    
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
fprintf('   PRESTRESSED BEAM ANALYSIS INPUT\n');
fprintf('========================================\n');

fprintf('\nBEAM GEOMETRY:\n');
fprintf('  Span length: %.1f in (%.2f ft)\n', beam.L, beam.L/12);
fprintf('  Analysis segments: %d\n', beam.num_segments);

fprintf('\nSECTION PROPERTIES:\n');
fprintf('  Area (A): %.2f in^2\n', section.A);
fprintf('  Moment of inertia (Ix): %.2f in^4\n', section.Ix);
fprintf('  Centroid from bottom (yc): %.2f in\n', section.yc);
fprintf('  Distance to top fiber (yt): %.2f in\n', section.yt);
fprintf('  Distance to bottom fiber (yb): %.2f in\n', section.yb);

fprintf('\nMATERIAL PROPERTIES:\n');
fprintf('  Concrete: f''c = %.1f ksi, Ec = %.0f ksi\n', materials.fc, materials.Ec);
fprintf('  Prestressing steel: fpu = %.0f ksi, Eps = %.0f ksi\n', materials.fpu, materials.Eps);
fprintf('  Mild steel: fy = %.0f ksi, Es = %.0f ksi\n', materials.fy, materials.Es);

fprintf('\nPRESTRESSING TENDONS:\n');
total_Aps = 0;
for i = 1:length(prestress.tendons)
    t = prestress.tendons{i};
    fprintf('  Tendon %d: Aps = %.2f in^2, fpi = %.1f ksi, Profile: %s\n', ...
        i, t.Aps, t.fpi, t.profile_type);
    fprintf('            Eccentricity: %.1f in (start) to %.1f in (mid) to %.1f in (end)\n', ...
        t.e_start, t.e_mid, t.e_end);
    fprintf('            Bonding: %s\n', t.bonding.type);
    total_Aps = total_Aps + t.Aps;
end
fprintf('  Total Aps: %.2f in^2\n', total_Aps);
fprintf('  Prestress losses: %.1f%%\n', prestress.losses * 100);

fprintf('\nNON-PRESTRESSED REINFORCEMENT:\n');
for i = 1:length(reinforcement.longitudinal)
    r = reinforcement.longitudinal{i};
    fprintf('  %s: y = %.1f in, As = %.2f in^2\n', ...
        r.type, r.y_position, sum(r.areas));
end

fprintf('\nAPPLIED LOADS:\n');
fprintf('  Distributed loads: %d segments\n', size(loads.distributed, 1));
fprintf('  Point loads: %d\n', size(loads.point, 1));
fprintf('  Support type: %s\n', loads.support_type);

fprintf('\n========================================\n');

end
