function [beam, section, materials, prestress, reinforcement, loads] = inputPrestressedBeam_Project1()
% INPUTPRESTRESSEDBEAM_PROJECT1 - Input data for Project 1 prestressed concrete beam
% Creates axial force, shear force, and bending moment diagrams
%
% Project 1 Specifications:
%   - Double-T (TT) section, symmetric about the y-axis
%   - Top flange: 120 in wide (from x=-60 to x=60), 2 in thick
%   - Two stems extending downward to y=-28 in
%   - Origin at the top surface of the top flange
%   - Prestressing tendons (to be defined per project specs)
%
% Section Vertices (x, y) [in] - counterclockwise from bottom-left:
%   (-60,  0)  -> Top-left of top flange
%   ( 60,  0)  -> Top-right of top flange
%   ( 60, -2)  -> Bottom-right of top flange
%   ( 33, -2)  -> Inner top-right of right stem
%   ( 32,-28)  -> Inner bottom-right of right stem
%   ( 28.25,-28) -> Outer bottom-right of right stem
%   ( 27.25,-2) -> Outer top-right of right stem (flange-stem junction)
%   (-27.25,-2) -> Outer top-left of left stem  (flange-stem junction)
%   (-28.25,-28) -> Outer bottom-left of left stem
%   (-32,-28)  -> Inner bottom-left of left stem
%   (-33, -2)  -> Inner top-left of left stem
%   (-60, -2)  -> Bottom-left of top flange
%
% Outputs:
%   beam         - Structure containing beam geometry and span information
%   section      - Structure containing cross-section geometry (polygonal)
%   materials    - Structure containing material properties
%   prestress    - Structure containing prestressing tendon information
%   reinforcement- Structure containing non-prestressed reinforcement
%   loads        - Structure containing applied loads

%% SECTION 1: BEAM GEOMETRY
beam.L            = 600;          % Total beam length (in) - UPDATE as needed
beam.num_segments = 100;          % Number of segments for analysis
beam.x = linspace(0, beam.L, beam.num_segments + 1);  % x-coordinates along beam

%% SECTION 2: CROSS-SECTION GEOMETRY
% Define vertices counterclockwise [x, y]
% Double-T section; origin at top surface of top flange
% y is POSITIVE upward (flange top = 0, stems go NEGATIVE downward)
%
%  _______________________________________________
% |                  TOP FLANGE                  |   y = 0  (top)
% |_______________________________________________|   y = -2 (flange soffit)
%        |           |           |           |
%     left gap    left stem   right stem   right gap
%        |           |           |           |
%        |___________|           |___________|   y = -28 (stem soffit)
%
section.vertices = [
    -60,      0;      % 1  Top-left of top flange
     60,      0;      % 2  Top-right of top flange
     60,     -2;      % 3  Bottom-right of top flange
     33,     -2;      % 4  Inner-top of right stem (flange edge)
     32,    -28;      % 5  Inner-bottom of right stem
     28.25, -28;      % 6  Outer-bottom of right stem
     27.25,  -2;      % 7  Outer-top of right stem (flange-stem junction)
    -27.25,  -2;      % 8  Outer-top of left stem  (flange-stem junction)
    -28.25, -28;      % 9  Outer-bottom of left stem
    -32,    -28;      % 10 Inner-bottom of left stem
    -33,     -2;      % 11 Inner-top of left stem (flange edge)
    -60,     -2;      % 12 Bottom-left of top flange
];

% Section properties (auto-calculated from vertices if left empty)
section.A   = [];   % Cross-sectional area (in^2)
section.Ix  = [];   % Moment of inertia about centroidal x-axis (in^4)
section.Iy  = [];   % Moment of inertia about centroidal y-axis (in^4)
section.yc  = [];   % Centroid y-coordinate from BOTTOM of section (in)
section.xc  = [];   % Centroid x-coordinate (in)  - should be 0 by symmetry
section.yt  = [];   % Distance from centroid to TOP fiber (in)
section.yb  = [];   % Distance from centroid to BOTTOM fiber (in)

%% SECTION 3: MATERIAL PROPERTIES
% --- Concrete ---
materials.fc  = 6.0;             % Compressive strength f'c (ksi)
materials.Ec  = 4700;            % Elastic modulus Ec (ksi)  -- UPDATE if needed
% Alternative: Ec = 57*sqrt(fc_psi) = 57*sqrt(6000) ~ 4415 ksi (ACI 318)
materials.fci = 4.8;             % Concrete strength at transfer f'ci (ksi)
%materials.fr  = 7.5 * sqrt(materials.fc * 1000) / 1000;  % Modulus of rupture (ksi)
materials.fr  = 3.0 * sqrt(materials.fc * 1000) / 1000;  % Modulus of rupture (ksi)

% --- Prestressing Steel ---
materials.fpu = 270;             % Ultimate tensile strength fpu (ksi)
materials.fpy = 0.9 * materials.fpu;  % Yield strength fpy (ksi)
materials.Eps = 28500;           % Elastic modulus Eps (ksi)

% --- Mild (Non-prestressed) Steel ---
materials.fy  = 60;              % Yield strength fy (ksi)
materials.Es  = 29000;           % Elastic modulus Es (ksi)

%% SECTION 4: PRESTRESSING TENDONS
% -----------------------------------------------------------------------
% NOTE: Update Aps, fpi, y_position (from BOTTOM of section), and
%       profile_type to match your Project 1 specifications.
% -----------------------------------------------------------------------

% Effective prestress after losses
% Assume initial jacking stress = 0.75 * fpu
fpi_init = 0.75 * materials.fpu;    % Initial prestress (ksi) = 202.5 ksi

% --- TENDON 1: Left stem ---
% Place at the centroid of the left stem bottom region (update y_position)
tendon1.Aps          = 2 * 0.153;   % Example: 2 x 0.5" dia strands = 0.306 in^2  -- UPDATE
tendon1.fpi          = fpi_init;    % Initial prestress (ksi)
tendon1.profile_type = 'linear'; % 'linear' | 'parabolic' | 'harped' | 'custom'
tendon1.x_position   = -29.875;     % x-centroid of left stem (in)
tendon1.y_position   = -22;         % y-coordinate from origin (in)
tendon1.e_start      = [];          % Will be computed from y_position & yc
tendon1.e_mid        = [];          % Will be computed from y_position & yc
tendon1.e_end        = [];          % Will be computed from y_position & yc
tendon1.drape_point  = [];          % Fraction of L (for harped profile only)
% Bonding
tendon1.bonding.type         = 'full';
tendon1.bonding.bonded_zones = [0, 1.0];  % Fully bonded

% --- TENDON 2: Right stem (mirror of Tendon 1) ---
tendon2.Aps          = 2 * 0.153;   % Same as Tendon 1  -- UPDATE
tendon2.fpi          = fpi_init;
tendon2.profile_type = 'linear';
tendon2.x_position   =  29.875;     % x-centroid of right stem (in)
tendon2.y_position   = -22;         % Same depth as Tendon 1
tendon2.e_start      = [];
tendon2.e_mid        = [];
tendon2.e_end        = [];
tendon2.drape_point  = [];
tendon2.bonding.type         = 'full';
tendon2.bonding.bonded_zones = [0, 1.0];

% Combine tendons (add or remove tendons as required by the project)
prestress.tendons = {tendon1, tendon2};
prestress.losses  = 0.15;           % Total prestress losses (fraction) -- 15% typical

%% SECTION 5: NON-PRESTRESSED REINFORCEMENT
% Update as required; leave empty if none specified in Project 1.

% Top reinforcement (compression steel)
rebar_top.y_position = [];
rebar_top.x_positions = [];
rebar_top.areas       = [];
rebar_top.type        = 'compression';

% Bottom reinforcement (tension steel)
rebar_bot.y_position  = [];
rebar_bot.x_positions = [];
rebar_bot.areas       = [];
rebar_bot.type        = 'tension';

% Shear reinforcement (stirrups)
stirrups.Av      = [];
stirrups.spacing = [];
stirrups.fy      = [];

% Combine reinforcement
reinforcement.longitudinal = {rebar_top, rebar_bot};
reinforcement.stirrups     = stirrups;

%% SECTION 6: APPLIED LOADS
% Self-weight (auto-calculated from section area and density if left empty)
loads.self_weight      = [];          % (kip/in) -- leave empty for auto-calc
loads.concrete_density = 150 / 1728; % kip/in^3  (150 pcf)

% Distributed loads: each row = [x_start, x_end, w_start, w_end] (kip/in)
% Negative w = downward.  Update or add rows as required.
loads.distributed = [
    % 0,    beam.L,  -0.10,  -0.10;   % Example: uniform superimposed dead load
];

% Point loads: each row = [x_position, P (kip), M (kip-in)]
% Positive P = upward; positive M = counterclockwise.
loads.point = [
    % beam.L/2,  -50,  0;   % Example: midspan concentrated load
];

% Support conditions
loads.support_type = 'simple';
loads.supports     = [0, beam.L];  % [x_left, x_right]

%% CALCULATE SECTION PROPERTIES
section = calculateSectionProperties(section);

%% PROCESS TENDON PROFILES
prestress = processTendonProfiles(prestress, beam, section);

%% DISPLAY SUMMARY
displayInputSummary(beam, section, materials, prestress, reinforcement, loads);

end

%% =========================================================================
%%  HELPER FUNCTIONS  (identical to HW1 version)
%% =========================================================================

function section = calculateSectionProperties(section)
% Calculate section properties from polygon vertices using the shoelace formula.

v = section.vertices;
n = size(v, 1);
v = [v; v(1,:)];   % Close polygon

% --- Area ---
A = 0;
for i = 1:n
    A = A + (v(i,1)*v(i+1,2) - v(i+1,1)*v(i,2));
end
A = abs(A) / 2;

% --- Centroid ---
Cx = 0; Cy = 0;
for i = 1:n
    f  = v(i,1)*v(i+1,2) - v(i+1,1)*v(i,2);
    Cx = Cx + (v(i,1) + v(i+1,1)) * f;
    Cy = Cy + (v(i,2) + v(i+1,2)) * f;
end
Cx = Cx / (6*A);
Cy = Cy / (6*A);

% --- Moments of inertia about centroidal axes ---
Ix = 0; Iy = 0;
for i = 1:n
    xi  = v(i,1)   - Cx;
    yi  = v(i,2)   - Cy;
    xi1 = v(i+1,1) - Cx;
    yi1 = v(i+1,2) - Cy;
    f   = xi*yi1 - xi1*yi;
    Ix  = Ix + (yi^2 + yi*yi1 + yi1^2) * f;
    Iy  = Iy + (xi^2 + xi*xi1 + xi1^2) * f;
end
Ix = abs(Ix) / 12;
Iy = abs(Iy) / 12;

% --- Assign ---
if isempty(section.A),   section.A  = A;  end
if isempty(section.Ix),  section.Ix = Ix; end
if isempty(section.Iy),  section.Iy = Iy; end
if isempty(section.xc),  section.xc = Cx; end
if isempty(section.yc),  section.yc = Cy; end

y_max = max(section.vertices(:,2));
y_min = min(section.vertices(:,2));
if isempty(section.yt), section.yt = y_max - section.yc; end
if isempty(section.yb), section.yb = section.yc - y_min; end

end

% --------------------------------------------------------------------------

function prestress = processTendonProfiles(prestress, beam, section)
% Generate eccentricity profile e(x) for each tendon.
% Eccentricity e = yc - y_tendon  (positive when tendon is BELOW centroid)

x  = beam.x;
L  = beam.L;
yc = section.yc;

for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};

    if isfield(tendon, 'y_position') && ~isempty(tendon.y_position)
        % Straight tendon at fixed y_position
        e_constant      = yc - tendon.y_position;
        tendon.e_start  = e_constant;
        tendon.e_mid    = e_constant;
        tendon.e_end    = e_constant;
        e               = ones(size(x)) * e_constant;
        y_tendon        = ones(size(x)) * tendon.y_position;

    else
        e_s = tendon.e_start;
        e_m = tendon.e_mid;
        e_e = tendon.e_end;

        switch tendon.profile_type
            case 'linear'
                e = e_s + (e_e - e_s) * x / L;

            case 'parabolic'
                % Passes through (0,e_s), (L/2,e_m), (L,e_e)
                a = 2*(e_s - 2*e_m + e_e) / L^2;
                b = (e_e - e_s)/L - a*L;
                c = e_s;
                e = a*x.^2 + b*x + c;

            case 'harped'
                x_drape = tendon.drape_point * L;
                e = zeros(size(x));
                m1 = x <= x_drape;
                m2 = x >  x_drape;
                e(m1) = e_s + (e_m - e_s) * x(m1) / x_drape;
                e(m2) = e_m + (e_e - e_m) * (x(m2) - x_drape) / (L - x_drape);

            case 'custom'
                if isfield(tendon, 'e_profile')
                    e = interp1(tendon.x_profile, tendon.e_profile, x, 'linear');
                else
                    error('Custom profile requires tendon.e_profile and tendon.x_profile');
                end

            otherwise
                error('Unknown profile_type: %s', tendon.profile_type);
        end

        y_tendon = yc - e;
    end

    prestress.tendons{i}.e      = e;
    prestress.tendons{i}.y      = y_tendon;
    prestress.tendons{i}.bonded = generateBondingMask(tendon.bonding, x, L);
end

end

% --------------------------------------------------------------------------

function bonded = generateBondingMask(bonding, x, L)
% Return logical array: true where tendon is bonded to concrete.

bonded = false(size(x));
switch bonding.type
    case 'full'
        bonded(:) = true;
    case 'unbonded'
        bonded(:) = false;
    case 'partial'
        for i = 1:size(bonding.bonded_zones, 1)
            xs = bonding.bonded_zones(i,1) * L;
            xe = bonding.bonded_zones(i,2) * L;
            bonded = bonded | (x >= xs & x <= xe);
        end
end

end

% --------------------------------------------------------------------------

function displayInputSummary(beam, section, materials, prestress, reinforcement, loads)
% Print a formatted summary of all input parameters.

fprintf('\n============================================\n');
fprintf('   PROJECT 1 - PRESTRESSED DOUBLE-T BEAM\n');
fprintf('============================================\n');

fprintf('\nBEAM GEOMETRY:\n');
fprintf('  Span length        : %.1f in  (%.2f ft)\n', beam.L, beam.L/12);
fprintf('  Analysis segments  : %d\n', beam.num_segments);

fprintf('\nSECTION PROPERTIES (computed):\n');
fprintf('  Area (A)           : %.3f in^2\n',  section.A);
fprintf('  Ix (centroidal)    : %.2f  in^4\n', section.Ix);
fprintf('  Centroid y (yc)    : %.4f in  (from origin)\n', section.yc);
fprintf('  Centroid x (xc)    : %.4f in  (should be ~0)\n', section.xc);
fprintf('  yt (top fiber)     : %.4f in\n', section.yt);
fprintf('  yb (bottom fiber)  : %.4f in\n', section.yb);

fprintf('\nMATERIAL PROPERTIES:\n');
fprintf('  Concrete : f''c = %.1f ksi (%5.0f psi),  Ec = %.0f ksi\n', ...
    materials.fc, materials.fc*1000, materials.Ec);
fprintf('  PS steel : fpu = %.0f ksi,  Eps = %.0f ksi\n', ...
    materials.fpu, materials.Eps);
fprintf('  Mild steel: fy  = %.0f ksi,  Es  = %.0f ksi\n', ...
    materials.fy, materials.Es);

fprintf('\nPRESTRESSING TENDONS:\n');
total_Aps = 0; total_P = 0;
for i = 1:length(prestress.tendons)
    t = prestress.tendons{i};
    P_init = t.Aps * t.fpi;
    fprintf('  Tendon %d:\n', i);
    fprintf('    Aps     = %.4f in^2\n', t.Aps);
    fprintf('    fpi     = %.1f  ksi\n', t.fpi);
    fprintf('    P_init  = %.1f  kips\n', P_init);
    if isfield(t, 'y_position') && ~isempty(t.y_position)
        fprintf('    y_pos   = %.2f in  (from origin)\n', t.y_position);
        fprintf('    e       = %.4f in  (constant)\n', t.e_start);
    else
        fprintf('    Profile = %s\n', t.profile_type);
        fprintf('    e       = %.2f -> %.2f -> %.2f in  (start/mid/end)\n', ...
            t.e_start, t.e_mid, t.e_end);
    end
    fprintf('    Bonding = %s\n', t.bonding.type);
    total_Aps = total_Aps + t.Aps;
    total_P   = total_P   + P_init;
end
fprintf('  Total Aps           : %.4f in^2\n', total_Aps);
fprintf('  Total P (initial)   : %.1f  kips\n', total_P);
fprintf('  Prestress losses    : %.1f %%\n', prestress.losses*100);
fprintf('  Total P (effective) : %.1f  kips\n', total_P*(1-prestress.losses));

fprintf('\nNON-PRESTRESSED REINFORCEMENT:\n');
has_rebar = false;
for i = 1:length(reinforcement.longitudinal)
    r = reinforcement.longitudinal{i};
    if ~isempty(r.y_position) && ~isempty(r.areas)
        fprintf('  %s: y = %.1f in,  As = %.3f in^2\n', ...
            r.type, r.y_position, sum(r.areas));
        has_rebar = true;
    end
end
if ~has_rebar, fprintf('  None specified\n'); end

fprintf('\nAPPLIED LOADS:\n');
if ~isempty(loads.distributed)
    fprintf('  Distributed loads : %d segment(s)\n', size(loads.distributed,1));
else
    fprintf('  Distributed loads : None\n');
end
if ~isempty(loads.point)
    fprintf('  Point loads       : %d\n', size(loads.point,1));
else
    fprintf('  Point loads       : None\n');
end
fprintf('  Support type      : %s\n', loads.support_type);
fprintf('  Supports at x     : %.1f in and %.1f in\n', ...
    loads.supports(1), loads.supports(2));

fprintf('\n============================================\n');
end
