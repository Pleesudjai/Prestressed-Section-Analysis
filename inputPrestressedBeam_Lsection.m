function [beam, section, materials, prestress, reinforcement, loads] = inputPrestressedBeam_Lsection()
% INPUTPRESTRESSEDBEAM_LSECTION - Input data for L-shaped prestressed beam
% Example using complex polygonal section similar to inputData.m sample
%
% This demonstrates how to define arbitrary section geometry

%% SECTION 1: BEAM GEOMETRY
beam.L = 480;                   % Total beam length (in) = 40 ft
beam.num_segments = 100;
beam.x = linspace(0, beam.L, beam.num_segments + 1);

%% SECTION 2: CROSS-SECTION GEOMETRY (L-Section)
% Define vertices counterclockwise [x, y] - origin at bottom-left
% This is an L-shaped section similar to the inputData.m example
section.vertices = [
    20.000,  240.000;   % Top of vertical leg
    0.000,   240.000;   % Top-left corner
    0,       0;         % Bottom-left (origin)
    150,     0;         % Bottom-right of horizontal leg
    150,     20;        % Top-right of horizontal leg
    20,      20;        % Inside corner
];

% Leave empty for automatic calculation
section.A = [];
section.Ix = [];
section.Iy = [];
section.yc = [];
section.xc = [];
section.yt = [];
section.yb = [];

%% SECTION 3: MATERIAL PROPERTIES
materials.fc = 9;               % ksi (high strength for UHPC/prestressed)
materials.Ec = 57 * sqrt(materials.fc * 1000);
materials.fci = 6.5;
materials.fr = 7.5 * sqrt(materials.fc * 1000) / 1000;

materials.fpu = 270;
materials.fpy = 0.9 * materials.fpu;
materials.Eps = 28500;

materials.fy = 60;
materials.Es = 29000;

%% SECTION 4: PRESTRESSING TENDONS
% Tendon 1 - In the vertical leg (parabolic)
tendon1.Aps = 1.5;
tendon1.fpi = 0.75 * materials.fpu;
tendon1.profile_type = 'parabolic';
tendon1.e_start = 50;           % Below centroid
tendon1.e_mid = 100;            % Maximum eccentricity at midspan
tendon1.e_end = 50;
tendon1.drape_point = [];

tendon1.bonding.type = 'partial';
tendon1.bonding.bonded_zones = [0, 0.15; 0.35, 0.65; 0.85, 1.0];

% Tendon 2 - In the horizontal leg (harped profile)
tendon2.Aps = 2.0;
tendon2.fpi = 0.72 * materials.fpu;
tendon2.profile_type = 'harped';
tendon2.e_start = 20;           
tendon2.e_mid = 80;
tendon2.e_end = 20;
tendon2.drape_point = 0.4;      % Drape point at 40% of span

tendon2.bonding.type = 'full';
tendon2.bonding.bonded_zones = [0, 1.0];

prestress.tendons = {tendon1, tendon2};
prestress.losses = 0.18;        % 18% total losses

%% SECTION 5: NON-PRESTRESSED REINFORCEMENT
% Along vertical leg (similar to inputData line1)
rebar_vertical.y_position = [];  % Will be distributed along height
rebar_vertical.x_positions = [2.5];
rebar_vertical.areas = [1.56];   % #11 bars
rebar_vertical.type = 'vertical_leg';

% Along horizontal leg (similar to inputData line2)
rebar_horizontal.y_position = 2.5;
rebar_horizontal.x_positions = [2.5, 12.5, 22.5, 32.5, 42.5, 52.5, 62.5, ...
                                72.5, 82.5, 92.5, 102.5, 112.5, 122.5, 132.5, 147.5];
rebar_horizontal.areas = 1.56 * ones(1, 15);
rebar_horizontal.type = 'tension';

% Inner corner reinforcement
rebar_corner.y_position = 17.5;
rebar_corner.x_positions = [17.5];
rebar_corner.areas = [1.56];
rebar_corner.type = 'corner';

stirrups.Av = 0.40;             % #4 stirrups, 2 legs
stirrups.spacing = 8;
stirrups.fy = 60;

reinforcement.longitudinal = {rebar_horizontal, rebar_corner};
reinforcement.stirrups = stirrups;

%% SECTION 6: APPLIED LOADS
loads.self_weight = [];
loads.concrete_density = 150/1728;

% Distributed loads
loads.distributed = [
    0,      480,    -0.08,  -0.08;   % Dead load
    48,     432,    -0.05,  -0.05;   % Live load (partial span)
];

% Point loads
loads.point = [
    160,    -30,    0;      % Point load at 1/3 span
    320,    -30,    0;      % Point load at 2/3 span
];

loads.support_type = 'simple';
loads.supports = [0, 480];

%% Calculate section properties
section = calculateSectionProperties(section);

%% Process tendon profiles
prestress = processTendonProfiles(prestress, beam, section);

%% Display summary
displayInputSummary(beam, section, materials, prestress, reinforcement, loads);

end

%% HELPER FUNCTIONS (same as main input file)
function section = calculateSectionProperties(section)
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

Ix = 0; Iy = 0;
for i = 1:n
    xi = v(i,1) - Cx;
    yi = v(i,2) - Cy;
    xi1 = v(i+1,1) - Cx;
    yi1 = v(i+1,2) - Cy;
    
    factor = (xi * yi1 - xi1 * yi);
    Ix = Ix + (yi^2 + yi*yi1 + yi1^2) * factor;
    Iy = Iy + (xi^2 + xi*xi1 + xi1^2) * factor;
end
Ix = abs(Ix) / 12;
Iy = abs(Iy) / 12;

if isempty(section.A), section.A = A; end
if isempty(section.Ix), section.Ix = Ix; end
if isempty(section.Iy), section.Iy = Iy; end
if isempty(section.xc), section.xc = Cx; end
if isempty(section.yc), section.yc = Cy; end

y_max = max(section.vertices(:,2));
y_min = min(section.vertices(:,2));
if isempty(section.yt), section.yt = y_max - section.yc; end
if isempty(section.yb), section.yb = section.yc - y_min; end
end

function prestress = processTendonProfiles(prestress, beam, section)
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
            e = e_s + (e_e - e_s) * x / L;
        case 'parabolic'
            a = 2 * (e_s - 2*e_m + e_e) / L^2;
            b = (e_e - e_s) / L - a * L;
            c = e_s;
            e = a * x.^2 + b * x + c;
        case 'harped'
            x_drape = tendon.drape_point * L;
            e = zeros(size(x));
            mask1 = x <= x_drape;
            mask2 = x > x_drape;
            e(mask1) = e_s + (e_m - e_s) * x(mask1) / x_drape;
            e(mask2) = e_m + (e_e - e_m) * (x(mask2) - x_drape) / (L - x_drape);
        otherwise
            error('Unknown profile type');
    end
    
    prestress.tendons{i}.e = e;
    prestress.tendons{i}.y = yc - e;
    prestress.tendons{i}.bonded = generateBondingMask(tendon.bonding, x, L);
end
end

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
            x_end = bonding.bonded_zones(i, 2) * L;
            bonded = bonded | (x >= x_start & x <= x_end);
        end
end
end

function displayInputSummary(beam, section, materials, prestress, reinforcement, loads)
fprintf('\n========================================\n');
fprintf('   L-SECTION PRESTRESSED BEAM INPUT\n');
fprintf('========================================\n');

fprintf('\nBEAM: L = %.1f in (%.2f ft)\n', beam.L, beam.L/12);
fprintf('\nSECTION: A = %.1f in², Ix = %.0f in⁴, yc = %.2f in\n', ...
    section.A, section.Ix, section.yc);
fprintf('\nMATERIALS: fc = %.1f ksi, fpu = %.0f ksi\n', materials.fc, materials.fpu);
fprintf('\nTENDONS: %d total, Total Aps = %.2f in²\n', ...
    length(prestress.tendons), sum(cellfun(@(t) t.Aps, prestress.tendons)));
fprintf('\n========================================\n');
end
