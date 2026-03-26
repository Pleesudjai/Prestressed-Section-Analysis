function [beam, section, materials, prestress, reinforcement, loads] = inputData_Tsection()
% INPUTDATA_TSECTION  — Simple T-beam test case
%
%  Purpose: verify the T-section branch (a > hf) in ultimateDesign_Tsection.m
%           with both mild steel (tension + compression) present.
%
%  Section:  inverted-T precast beam
%    b   = 36 in  (effective flange width)
%    bw  = 12 in  (web width)
%    hf  =  4 in  (flange thickness)
%    h   = 36 in  (total depth)
%
%  Coordinate origin: y = 0 at BOTTOM of web, y = 36 at top of flange.
%
%  ASCII cross-section (y upward):
%
%  y=36  +---------------------------+  <- top of flange
%  y=32  +--------+         +--------+  <- flange soffit
%                 |  web    |
%  y= 4      ->   |  o  o   |  <- PS tendons at y=4 in
%  y= 0           +---------+  <- bottom of web
%
%  Vertices (counterclockwise, starting at bottom-left of web):
%   (-6, 0)  -> (6, 0)  -> (6, 32) -> (18, 32) -> (18, 36) ->
%   (-18, 36) -> (-18, 32) -> (-6, 32) -> back to start
%
%  Steel:
%    Aps  = 2.40 in^2  (single equivalent tendon, y = 4 in from bottom = dp=32 from top)
%    As   = 2.00 in^2  tension mild steel, y = 3 in from bottom (d = 33 in from top)
%    As'  = 0.50 in^2  compression mild steel, y = 33 in from bottom (d' = 3 in from top)
%
%  Span: 40 ft = 480 in, simply supported.

%% SECTION 1: BEAM GEOMETRY
beam.L            = 40 * 12;       % 480 in (40 ft)
beam.num_segments = 500;

beam.x_plot_fractions = [0, 0.25, 0.5, 1.0];
beam.x = linspace(0, beam.L, beam.num_segments + 1);
x_plot_exact = beam.x_plot_fractions * beam.L;
beam.x = unique([beam.x, x_plot_exact]);

%% SECTION 2: CROSS-SECTION GEOMETRY (T-section, y=0 at web bottom)
%  Vertices counterclockwise:
section.vertices = [
    -6,   0;    % 1  bottom-left of web
     6,   0;    % 2  bottom-right of web
     6,  32;    % 3  top-right of web  (flange soffit)
    18,  32;    % 4  outer-right of flange soffit
    18,  36;    % 5  top-right of flange
   -18,  36;    % 6  top-left of flange
   -18,  32;    % 7  outer-left of flange soffit
    -6,  32;    % 8  top-left of web  (flange soffit)
];

% Leave empty — auto-calculated by shoelace in ultimateDesign_Tsection.m
section.A   = [];
section.Ix  = [];
section.Iy  = [];
section.yc  = [];
section.xc  = [];
section.yt  = [];
section.yb  = [];

%% SECTION 3: MATERIAL PROPERTIES
materials.fc  = 4.0;             % f'c  (ksi)
materials.Ec  = 57 * sqrt(4000) / 1000;   % ~3.605 ksi  (ACI 318 Ec = 57√f'c_psi)
materials.fci = 3.2;             % f'ci (ksi) — 80% of f'c (typical at release)
materials.fr  = 7.5 * sqrt(materials.fc * 1000) / 1000;   % modulus of rupture (ksi)
materials.code_edition = 'ACI-318-19';

% Allowable stresses at TRANSFER
materials.f_ci_allow     =  0.60 * materials.fci;
materials.f_ti_allow     = -3.0  * sqrt(materials.fci * 1000) / 1000;
materials.f_ti_allow_end = -6.0  * sqrt(materials.fci * 1000) / 1000;
materials.f_ci_allow_end =  0.70 * materials.fci;   % ACI-318-19

% Allowable stresses at SERVICE
materials.f_cs_allow_sust  =  0.45 * materials.fc;
materials.f_cs_allow_total =  0.60 * materials.fc;
materials.f_ts_allow       = -12.0 * sqrt(materials.fc * 1000) / 1000;   % Class C
materials.f_tu_allow       =  -7.5 * sqrt(materials.fc * 1000) / 1000;   % Class U (ACI-318-19)

% Prestressing steel  (low-relaxation — fpy/fpu = 0.90 → gamma_p = 0.28)
materials.fpu = 270;
materials.fpy = 0.90 * materials.fpu;   % = 243 ksi
materials.Eps = 28500;

% Mild steel
materials.fy  = 60;
materials.Es  = 29000;

%% SECTION 4: PRESTRESSING TENDONS
%  Single equivalent tendon with Aps = 2.40 in^2 at y = 4 in from bottom
%    fpi  = 0.70 * fpu = 0.70 * 270 = 189 ksi
%    Pi   = Aps * fpi = 2.40 * 189 = 453.6 kip
%    fse  = 0.85 * 189 = 160.65 ksi  (after 15% losses)
fpi_val = 0.70 * materials.fpu;   % 189 ksi initial stress

tendon1.Aps          = 2.40;         % in^2  (equivalent to 16 × 0.150 or similar)
tendon1.fpi          = fpi_val;      % 189 ksi
tendon1.profile_type = 'linear';     % straight tendon
tendon1.x_position   = 0;           % centred in web (cosmetic only)
tendon1.y_position   = 4;           % y = 4 in from bottom  →  dp = 36-4 = 32 in from top
tendon1.e_start      = [];
tendon1.e_mid        = [];
tendon1.e_end        = [];
tendon1.drape_point  = [];
tendon1.bonding.type         = 'full';
tendon1.bonding.bonded_zones = [0, 1.0];

prestress.tendons = {tendon1};
prestress.losses  = 0.15;            % 15% total losses

%% SECTION 5: NON-PRESTRESSED REINFORCEMENT
%  TENSION MILD STEEL  (As = 2.00 in^2 near bottom)
%    y_position = 3 in from bottom   →   d = 36 - 3 = 33 in from top
rebar_bot.y_position  = 3;           % in from bottom (y-origin)
rebar_bot.x_positions = [beam.L/2];  % midspan (representative; used for display only)
rebar_bot.areas       = [2.00];      % in^2 total
rebar_bot.type        = 'tension';

%  COMPRESSION MILD STEEL  (As' = 0.50 in^2 near top)
%    y_position = 33 in from bottom  →   d' = 36 - 33 = 3 in from top
rebar_top.y_position  = 33;          % in from bottom
rebar_top.x_positions = [beam.L/2];
rebar_top.areas       = [0.50];      % in^2 total
rebar_top.type        = 'compression';

stirrups.Av      = [];
stirrups.spacing = [];
stirrups.fy      = [];

reinforcement.longitudinal = {rebar_top, rebar_bot};
reinforcement.stirrups     = stirrups;

%% SECTION 6: APPLIED LOADS
loads.self_weight      = [];
loads.concrete_density = 150 / 1728;   % lb/in^3

% SDL: 25 psf superimposed dead (e.g. floor finish)
%   w_SDL = 25 psf * 3 ft tributary / 1000 = 0.075 kip/ft = 0.00625 kip/in
w_SDL = 25 * 3 / (12 * 1000);    % kip/in  (25 psf * 3-ft tributary)

% LL: 100 psf live load, 3 ft tributary
%   w_LL = 100 psf * 3 ft / 1000 = 0.300 kip/ft = 0.025 kip/in
w_LL  = 100 * 3 / (12 * 1000);   % kip/in

loads.distributed = [
    0, beam.L, -w_SDL, -w_SDL;
];
loads.distributed_live = [
    0, beam.L, -w_LL, -w_LL;
];
loads.point = [];

loads.support_type = 'simple';
loads.supports     = [0, beam.L];

%% CALCULATE SECTION PROPERTIES
section = calculateSectionProperties(section);

%% PROCESS TENDON PROFILES
prestress = processTendonProfiles(prestress, beam, section);

%% DISPLAY SUMMARY
displayInputSummary(beam, section, materials, prestress, reinforcement, loads);

end

%% =========================================================================
%%  HELPER FUNCTIONS  (identical to inputData.m)
%% =========================================================================

function section = calculateSectionProperties(section)
v = section.vertices;
n = size(v, 1);
v = [v; v(1,:)];
A = 0;
for i = 1:n
    A = A + (v(i,1)*v(i+1,2) - v(i+1,1)*v(i,2));
end
A = abs(A) / 2;
Cx = 0; Cy = 0;
for i = 1:n
    f  = v(i,1)*v(i+1,2) - v(i+1,1)*v(i,2);
    Cx = Cx + (v(i,1) + v(i+1,1)) * f;
    Cy = Cy + (v(i,2) + v(i+1,2)) * f;
end
Cx = Cx / (6*A);
Cy = Cy / (6*A);
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

function prestress = processTendonProfiles(prestress, beam, section)
x  = beam.x;
L  = beam.L;
yc = section.yc;
for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};
    if isfield(tendon, 'y_position') && ~isempty(tendon.y_position)
        e_constant      = yc - tendon.y_position;
        tendon.e_start  = e_constant;
        tendon.e_mid    = e_constant;
        tendon.e_end    = e_constant;
        e               = ones(size(x)) * e_constant;
        y_tendon        = ones(size(x)) * tendon.y_position;
    else
        if isfield(tendon, 'y_start') && ~isempty(tendon.y_start)
            tendon.e_start = yc - tendon.y_start;
        end
        if isfield(tendon, 'y_mid') && ~isempty(tendon.y_mid)
            tendon.e_mid = yc - tendon.y_mid;
        end
        if isfield(tendon, 'y_end') && ~isempty(tendon.y_end)
            tendon.e_end = yc - tendon.y_end;
        end
        e_s = tendon.e_start;
        e_m = tendon.e_mid;
        e_e = tendon.e_end;
        switch tendon.profile_type
            case 'linear'
                e = e_s + (e_e - e_s) * x / L;
            case 'parabolic'
                a = 2*(e_s - 2*e_m + e_e) / L^2;
                b_c = (e_e - e_s)/L - a*L;
                c_c = e_s;
                e = a*x.^2 + b_c*x + c_c;
            case 'harped'
                x_drape = tendon.drape_point * L;
                e = zeros(size(x));
                m1 = x <= x_drape;
                m2 = x >  x_drape;
                e(m1) = e_s + (e_m - e_s) * x(m1) / x_drape;
                e(m2) = e_m + (e_e - e_m) * (x(m2) - x_drape) / (L - x_drape);
            case 'custom'
                if isfield(tendon, 'y_profile') && ~isempty(tendon.y_profile)
                    e_profile = yc - tendon.y_profile;
                    e = interp1(tendon.x_profile, e_profile, x, 'linear');
                elseif isfield(tendon, 'e_profile')
                    e = interp1(tendon.x_profile, tendon.e_profile, x, 'linear');
                else
                    error('Custom profile requires y_profile or e_profile');
                end
            otherwise
                error('Unknown profile_type: %s', tendon.profile_type);
        end
        y_tendon = yc - e;
    end
    prestress.tendons{i}.e       = e;
    prestress.tendons{i}.y       = y_tendon;
    prestress.tendons{i}.e_start = e(1);
    prestress.tendons{i}.e_mid   = e(round(length(e)/2));
    prestress.tendons{i}.e_end   = e(end);
    prestress.tendons{i}.bonded  = generateBondingMask(tendon.bonding, x, L);
end
end

function bonded = generateBondingMask(bonding, x, L)
bonded = false(size(x));
switch bonding.type
    case 'full',    bonded(:) = true;
    case 'unbonded', bonded(:) = false;
    case 'partial'
        for i = 1:size(bonding.bonded_zones, 1)
            xs = bonding.bonded_zones(i,1) * L;
            xe = bonding.bonded_zones(i,2) * L;
            bonded = bonded | (x >= xs & x <= xe);
        end
end
end

function displayInputSummary(beam, section, materials, prestress, reinforcement, loads)
fprintf('\n============================================\n');
fprintf('   T-SECTION TEST CASE — VERIFY T-SECTION BRANCH\n');
fprintf('   b=36 in  bw=12 in  hf=4 in  h=36 in\n');
fprintf('============================================\n');
fprintf('\nBEAM GEOMETRY:\n');
fprintf('  Span length        : %.1f in  (%.2f ft)\n', beam.L, beam.L/12);
fprintf('  Analysis segments  : %d\n', beam.num_segments);
fprintf('\nSECTION PROPERTIES (computed):\n');
fprintf('  Area (A)           : %.3f in^2\n',  section.A);
fprintf('  Ix (centroidal)    : %.2f  in^4\n', section.Ix);
fprintf('  Centroid yc        : %.4f in  (from bottom)\n', section.yc);
fprintf('  yt (centroid->top) : %.4f in\n', section.yt);
fprintf('  yb (centroid->bot) : %.4f in\n', section.yb);
fprintf('\nMATERIAL PROPERTIES:\n');
fprintf('  f''c = %.1f ksi,  fpu = %.0f ksi,  fy = %.0f ksi\n', ...
    materials.fc, materials.fpu, materials.fy);
fprintf('\nPRESTRESSING TENDONS:\n');
for i = 1:length(prestress.tendons)
    t = prestress.tendons{i};
    fprintf('  Tendon %d: Aps=%.3f in^2,  fpi=%.1f ksi,  y=%.1f in from bottom\n', ...
        i, t.Aps, t.fpi, t.y_position);
end
total_Aps = sum(cellfun(@(t) t.Aps, prestress.tendons));
total_Pi  = sum(cellfun(@(t) t.Aps*t.fpi, prestress.tendons));
fprintf('  Total Aps = %.4f in^2,  Total Pi = %.1f kip\n', total_Aps, total_Pi);
fprintf('  Losses = %.0f%%,  Total Pe = %.1f kip\n', prestress.losses*100, total_Pi*(1-prestress.losses));
fprintf('\nNON-PRESTRESSED REINFORCEMENT:\n');
for i = 1:length(reinforcement.longitudinal)
    r = reinforcement.longitudinal{i};
    if ~isempty(r.y_position) && ~isempty(r.areas)
        fprintf('  %s: y = %.1f in from bottom,  A = %.3f in^2\n', ...
            r.type, r.y_position, sum(r.areas));
    end
end
fprintf('\n============================================\n');
end
