%% ULTIMATE_DESIGN  — Nominal Moment Capacity (ACI 318-19 / CEE530)
%  Project 1: Double-T Prestressed Beam
%
%  HOW TO RUN:
%    1. Set MATLAB working directory to project_Project2/
%    2. Run this script (F5 or >> ultimateDesign)
%    3. Figures + report saved to  output/
%
%  SIGN CONV.:   compression = positive, tension = negative
%
%  KEY FORMULAS:
%    fps = fpu·[1 − (γp/β1)·ρp·fpu/f'c]    (bonded, prestress only)
%    a_trial = Aps·fps / (0.85·f'c·b)
%    If a ≤ hf → rectangular;  else → T-section (recompute with bw)
%    φMn ≥ Mu  AND  φMn ≥ 1.2·Mcr
%    εt = (dp−c)/c · 0.003 ≥ 0.004 (ACI §9.3.3)

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));  % shared helpers
clear; clc; close all;

%% ======================================================================
%%  USER TOGGLE — Design code edition
%% ======================================================================
%  'CEE530'     : CEE 530 course parameters
%                   fpy/fpu = 0.85  →  γp = 0.40  (stress-relieved strand)
%                   Matches Group 3 PDF submission
%  'ACI-318-19' : Current ACI standard
%                   fpy/fpu = 0.90  →  γp = 0.28  (low-relaxation strand)
%                   Matches inputData.m  (fpy = 0.9·fpu = 243 ksi)
design_code = 'ACI-318-19';   % T-section test: low-relax strand -> gamma_p = 0.28

%% ======================================================================
%%  USER TOGGLE — Number of strands (design Aps override)
%% ======================================================================
%  inputData.m defines 4 strands (original project).
%  Ultimate design requires 13 strands to satisfy φMn ≥ Mu.
%  Set n_strands_design = [] to use inputData.m value (4 strands).
n_strands_design = [];    % T-section test: use Aps directly from inputData_Tsection
Aps_per_strand   = 0.153; % in^2 per 1/2-in strand

output_dir = 'output';
if ~exist(output_dir, 'dir');  mkdir(output_dir);  end

%% ======================================================================
%%  STEP 0 — Load input data
%% ======================================================================
fprintf('Loading input data...\n');
[beam, section, materials, prestress, reinforcement, loads] = inputData_Tsection();

%% ======================================================================
%%  STEP 1 — Compute section properties from vertices (shoelace)
%% ======================================================================
vx = section.vertices(:,1);
vy = section.vertices(:,2);
n_v = length(vx);

% Shoelace formula
A_sec = 0;
Qx    = 0;
for k = 1:n_v
    kn = mod(k, n_v) + 1;
    cross_k = vx(k)*vy(kn) - vx(kn)*vy(k);
    A_sec = A_sec + cross_k;
    Qx    = Qx    + (vy(k) + vy(kn)) * cross_k;
end
A_sec = abs(A_sec) / 2;
yc    = abs(Qx) / (6 * A_sec);   % centroid from bottom

% Moment of inertia (about centroidal axis)
Ic = 0;
for k = 1:n_v
    kn = mod(k, n_v) + 1;
    cross_k = vx(k)*vy(kn) - vx(kn)*vy(k);
    Ic = Ic + (vy(k)^2 + vy(k)*vy(kn) + vy(kn)^2) * cross_k;
end
Ic  = abs(Ic) / 12 - A_sec * yc^2;
h   = max(vy);       % total height
yt  = h - yc;        % distance from centroid to TOP fiber
yb  = yc;            % distance from centroid to BOTTOM fiber
St  = Ic / yt;       % top section modulus
Sb  = Ic / yb;       % bottom section modulus

fprintf('\n--- SECTION PROPERTIES ---\n');
fprintf('  A  = %.3f in^2\n', A_sec);
fprintf('  yc = %.4f in  (from bottom)\n', yc);
fprintf('  yt = %.4f in  (centroid to top)\n', yt);
fprintf('  yb = %.4f in  (centroid to bottom)\n', yb);
fprintf('  Ic = %.1f in^4\n', Ic);
fprintf('  St = %.1f in^3\n', St);
fprintf('  Sb = %.1f in^3\n', Sb);
fprintf('  h  = %.1f in\n', h);

%% ======================================================================
%%  STEP 2 — T-section geometry (effective flange + web)
%% ======================================================================
% For the Double-T:
%   b  = full flange width = 120 in
%   hf = flange thickness  = 2 in   (from y=26 to y=28)
%   bw = total stem width  = 2 × (28.25−27.25) + 2 × (32−28.25) = 2×1+2×3.75 = no
%      Exact stem widths from vertices:
%        Left stem outer: x from -28.25 to -27.25  (1 in at top)
%                         x from -32    to -28.25  (3.75 in at bottom)
%        → average/effective = 32-27.25 = 4.75 in per stem (use rectangular approx.)
%   For ACI fps / a calculations, use the RECTANGULAR stem width.
%   The stem sides are slightly tapered; use the AVERAGE stem width:
%      Each stem average width ≈ (outerTop−innerTop)/2 + ... but simpler:
%      Inner face at y=0: x = ±32,   ±28.25
%      Inner face at y=26: x = ±33,  ±27.25
%      Average inner width of each stem = (32+33)/2 - (28.25+27.25)/2 = 32.5 - 27.75 = 4.75 in
%   bw = 2 × 4.75 = 9.50 in  (total, both stems)
%
%  NOTE: Use the wider average so a conservative (smaller) fps results.
b  = 36;         % T-section test: flange width (in)
hf = 4;          % T-section test: flange thickness (in)  [y = 32 to 36]
bw = 12;         % T-section test: web width (in)

fprintf('\n--- T-SECTION GEOMETRY ---\n');
fprintf('  b  (flange) = %.1f in\n', b);
fprintf('  hf (flange) = %.1f in\n', hf);
fprintf('  bw (webs)   = %.2f in  (2 stems × 4.75 in avg)\n', bw);

%% ======================================================================
%%  STEP 3 — Material factors
%% ======================================================================
fc     = materials.fc;     % f'c  (ksi)
fpu    = materials.fpu;    % 270 ksi
losses = prestress.losses; % 0.15
eta    = 1 - losses;       % 0.85

% ACI β1 (Table 22.2.2.4.3)
beta1 = 0.85 - 0.05 * (fc - 4);
beta1 = max(0.65, min(0.85, beta1));   % ACI bounds: [0.65, 0.85]

% γp — controlled by design_code toggle
if strcmp(design_code, 'CEE530')
    % CEE530 course: stress-relieved strand, fpy/fpu = 0.85
    fpy     = 0.85 * fpu;   % = 229.5 ksi
    gamma_p = 0.40;
    code_label = 'CEE530 (stress-relieved, fpy/fpu=0.85, γp=0.40)';
else
    % ACI 318-19: low-relaxation strand, fpy/fpu = 0.90
    fpy     = materials.fpy;   % = 243 ksi (from inputData)
    gamma_p = 0.28;
    code_label = 'ACI-318-19 (low-relaxation, fpy/fpu=0.90, γp=0.28)';
end

fprintf('\n--- MATERIAL FACTORS  [%s] ---\n', code_label);
fprintf('  f''c     = %.1f ksi\n', fc);
fprintf('  fpu     = %.0f ksi\n', fpu);
fprintf('  fpy     = %.1f ksi  (fpy/fpu = %.2f)\n', fpy, fpy/fpu);
fprintf('  eta     = %.2f  (1 - %.0f%% losses)\n', eta, losses*100);
fprintf('  beta1   = %.2f  (ACI Table 22.2.2.4.3)\n', beta1);
fprintf('  gamma_p = %.2f  (ACI Table 22.3.2.1)\n', gamma_p);

%% ======================================================================
%%  STEP 4 — Tendon effective prestress check  (fse ≥ 0.5·fpu)
%% ======================================================================
% fse = effective steel stress after losses
%     = fpi × eta   (simplified — initial prestress × loss factor)
tendons = prestress.tendons;
n_tend  = numel(tendons);

% Total Aps and force-weighted centroid at midspan
Aps_total = 0;
Fps_total = 0;   % sum of Aps_i × fse_i (for centroid calc)
y_weighted = 0;
fse_check_pass = true;

for k = 1:n_tend
    td     = tendons{k};
    Aps_k  = td.Aps;
    fpi_k  = td.fpi;
    fse_k  = fpi_k * eta;   % effective prestress (ksi)

    % y of tendon k at midspan
    x_mid = beam.L / 2;
    switch td.profile_type
        case 'linear'
            y_k = td.y_position;
        case {'harped','parabolic'}
            y_k = td.y_position;   % assume constant if linear
        case 'custom'
            y_k = interp1(td.x_profile, td.y_profile, x_mid, 'linear');
        otherwise
            y_k = td.y_position;
    end

    Aps_total  = Aps_total  + Aps_k;
    Fps_total  = Fps_total  + Aps_k * fse_k;
    y_weighted = y_weighted + Aps_k * y_k;

    if fse_k < 0.5 * fpu
        fse_check_pass = false;
        fprintf('  WARNING: Tendon %d  fse=%.1f ksi < 0.5*fpu=%.1f ksi\n', ...
                k, fse_k, 0.5*fpu);
    end
end

y_tendon_mid = y_weighted / Aps_total;   % force-weighted centroid at midspan (y from bottom)
fse_avg      = Fps_total  / Aps_total;   % average effective prestress (ksi)

% --- Aps override: replace strand count if n_strands_design is set ---
if ~isempty(n_strands_design)
    Aps_from_input = Aps_total;
    Aps_total = n_strands_design * Aps_per_strand;
    fprintf('\n  [OVERRIDE] Aps: inputData has %.3f in^2 (%d strands)\n', ...
            Aps_from_input, round(Aps_from_input / Aps_per_strand));
    fprintf('             Using design value: %d strands × %.3f in^2 = %.3f in^2\n', ...
            n_strands_design, Aps_per_strand, Aps_total);
    % y_tendon_mid unchanged — all strands at same dp
end

Fe = Aps_total * fse_avg;   % total effective prestress force (kips)

if fse_check_pass
    fprintf('\n  fse = %.1f ksi  ≥  0.5·fpu = %.1f ksi  → ACI fps formula is valid ✓\n', ...
            fse_avg, 0.5*fpu);
end

%% ======================================================================
%%  STEP 5 — Effective depth dp at midspan (correct formula)
%%
%%  dp = distance from EXTREME COMPRESSION FIBER to centroid of tension steel
%%     = h − y_tendon          (y measured from bottom)
%%  Equivalently:
%%     dp = yt + e_mid          where e = yc − y_tendon (eccentricity)
%% ======================================================================
dp  = h - y_tendon_mid;        % = yt + (yc - y_tendon_mid)  [in]
e_mid = yc - y_tendon_mid;     % eccentricity at midspan [in] (positive = below centroid)

fprintf('\n--- EFFECTIVE DEPTH dp (at midspan) ---\n');
fprintf('  h              = %.4f in\n', h);
fprintf('  y_tendon_mid   = %.4f in  (force-weighted centroid, from bottom)\n', y_tendon_mid);
fprintf('  yc             = %.4f in\n', yc);
fprintf('  e_mid          = yc - y_tendon = %.4f in\n', e_mid);
fprintf('  dp = h - y_tendon = %.4f in\n', dp);

%% ======================================================================
%%  STEP 6 — fps (bonded tendons, ACI 318-19 Eq. 22.3.2.1 full form)
%%
%%  Full form (with mild steel tension + compression):
%%    fps = fpu · [1 − (γp/β1) · (ρp·fpu/f'c + (d/dp)·ω − (d'/dp)·ω')]
%%  where  ω  = ρs · fy/f'c  (tension mild steel)
%%         ω' = ρs'· fy/f'c  (compression mild steel)
%%
%%  When As = 0 and As' = 0:  ω = ω' = 0  → reduces to prestress-only form.
%%  ρp computed with b for flange check, then bw if T-section needed.
%% ======================================================================

% --- Extract mild steel from reinforcement struct ---
fy   = materials.fy;   % ksi
As   = 0;   d_s  = 0;   % tension mild steel (sum over bars)
Asc  = 0;   d_sc = 0;   % compression mild steel (sum over bars)
for ii = 1:numel(reinforcement.longitudinal)
    r = reinforcement.longitudinal{ii};
    if isempty(r.areas), continue; end
    Ar = sum(r.areas);
    if Ar == 0, continue; end
    if strcmp(r.type, 'tension')
        As  = As  + Ar;
        d_s = d_s + Ar * (h - r.y_position);   % depth from top
    else
        Asc  = Asc  + Ar;
        d_sc = d_sc + Ar * (h - r.y_position);  % depth from top
    end
end
if As  > 0, d_s  = d_s  / As;  end   % weighted depth
if Asc > 0, d_sc = d_sc / Asc; end

% omega and omega' (force ratios — zero when no mild steel)
omega   = (As  / (b * dp)) * fy / fc;   % uses b (flange) for trial
omega_c = (Asc / (b * dp)) * fy / fc;

% Pre-initialise T-section variables (used in report even if rectangular branch taken)
rho_pw = NaN;  A_fl = NaN;  A_web = NaN;  ybar = NaN;

% --- Step 6a: trial fps using full flange width b ---
rho_p_trial = Aps_total / (b * dp);
fps_trial   = fpu * (1 - (gamma_p / beta1) * ...
              (rho_p_trial * fpu / fc  +  (d_s/dp)*omega  -  (d_sc/dp)*omega_c));

% --- Step 6b: trial compression block depth a_trial ---
%   Full equilibrium: Aps*fps + As*fy - Asc*fy = 0.85*f'c*b*a
a_trial = (Aps_total * fps_trial + As*fy - Asc*fy) / (0.85 * fc * b);

fprintf('\n--- fps CALCULATION ---\n');
fprintf('  [Trial — rectangular, b = %.0f in]\n', b);
fprintf('  rho_p     = Aps/(b·dp) = %.4f/(%.0f×%.4f) = %.6f\n', ...
        Aps_total, b, dp, rho_p_trial);
fprintf('  fps_trial = %.2f ksi  (see Part 6 derivation above)\n', fps_trial);
fprintf('  a_trial   = (Aps·fps + As·fy - As''·fy) / (0.85·f''c·b)\n');
fprintf('            = (%.4f×%.2f + %.4f×%.0f - %.4f×%.0f) / (0.85×%.1f×%.0f)\n', ...
        Aps_total, fps_trial, As, fy, Asc, fy, fc, b);
fprintf('            = %.4f in\n', a_trial);

if a_trial <= hf
    %----------------------------------------------------------------------
    % RECTANGULAR SECTION  (compression block lies entirely in flange)
    %----------------------------------------------------------------------
    a   = a_trial;
    fps = fps_trial;
    %  Mn = Aps·fps·(dp - a/2) + As·fy·(d - a/2) - Asc·fy·(d' - a/2)
    Mn  = Aps_total*fps*(dp - a/2) ...
        + As *fy *(d_s  - a/2) ...
        - Asc*fy *(d_sc - a/2);
    section_type = 'Rectangular (a <= hf)';

    fprintf('\n  a_trial = %.4f in  <=  hf = %.1f in  -->  RECTANGULAR section\n', a, hf);
    fprintf('  Compression block lies entirely within the flange. b_eff = %.0f in.\n', b);
    fprintf('\n--- NOMINAL MOMENT Mn ---\n');
    fprintf('  fps  = %.2f ksi\n', fps);
    fprintf('  a    = %.4f in\n', a);
    fprintf('  Formula: Mn = Aps*fps*(dp - a/2) + As*fy*(d - a/2) - As''*fy*(d'' - a/2)\n');
    fprintf('         [With As = %.4f in^2,  As'' = %.4f in^2]\n', As, Asc);
    fprintf('  Mn   = %.4f*%.2f*(%.4f - %.4f/2)\n', Aps_total, fps, dp, a);
    if As > 0
        fprintf('       + %.4f*%.0f*(%.4f - %.4f/2)\n', As, fy, d_s, a);
    end
    if Asc > 0
        fprintf('       - %.4f*%.0f*(%.4f - %.4f/2)\n', Asc, fy, d_sc, a);
    end
    fprintf('       = %.1f kip-in  (%.1f kip-ft)\n', Mn, Mn/12);

else
    %----------------------------------------------------------------------
    % T-SECTION  (neutral axis extends into web)
    %----------------------------------------------------------------------
    fprintf('\n  a_trial = %.4f in  >  hf = %.1f in  -->  T-SECTION required\n', a_trial, hf);
    fprintf('  Compression block penetrates into the web. Must split into flange + web.\n');

    % Recompute fps using web steel ratio (bw replaces b)
    rho_pw  = Aps_total / (bw * dp);
    omega_w = (As  / (bw * dp)) * fy / fc;
    omega_cw= (Asc / (bw * dp)) * fy / fc;
    fps    = fpu * (1 - (gamma_p / beta1) * ...
             (rho_pw * fpu / fc  +  (d_s/dp)*omega_w  -  (d_sc/dp)*omega_cw));

    % Flange overhang compression force (flanges only, not the web strip)
    Ff = 0.85 * fc * (b - bw) * hf;

    % Web compression block depth from equilibrium:
    %   Aps*fps + As*fy - Asc*fy = Ff + 0.85*f'c*bw*a
    a = (Aps_total*fps + As*fy - Asc*fy - Ff) / (0.85 * fc * bw);

    % Compression block geometry
    A_fl  = (b - bw) * hf;            % flange overhang area
    A_web = bw * a;                    % web compression area
    % Centroid of compression resultant from top fiber:
    ybar = (A_fl*(hf/2) + A_web*(a/2)) / (A_fl + A_web);

    % Mn = sum of tension forces × their lever arms to ybar
    Mn = Aps_total*fps*(dp  - ybar) ...
       + As *fy *(d_s  - ybar) ...
       - Asc*fy *(d_sc - ybar);
    section_type = 'T-Section (a > hf)';

    fprintf('\n--- T-SECTION RECOMPUTATION ---\n');
    fprintf('  Use web width bw = %.2f in for fps and web block depth.\n\n', bw);
    fprintf('  rho_pw = Aps/(bw*dp) = %.4f/(%.2f*%.4f) = %.6f\n', ...
            Aps_total, bw, dp, rho_pw);
    fprintf('  fps    = %.2f ksi  (full formula with bw)\n\n', fps);
    fprintf('  Compression zone — split into FLANGE overhangs + WEB:\n');
    fprintf('    Flange overhang force:  Ff = 0.85*f''c*(b-bw)*hf\n');
    fprintf('                               = 0.85*%.1f*(%.0f-%.2f)*%.1f = %.2f kip\n', ...
            fc, b, bw, hf, Ff);
    fprintf('    Web block depth (from equilibrium):\n');
    fprintf('      a = (Aps*fps + As*fy - As''*fy - Ff) / (0.85*f''c*bw)\n');
    fprintf('        = (%.4f*%.2f + %.4f*%.0f - %.4f*%.0f - %.2f) / (0.85*%.1f*%.2f)\n', ...
            Aps_total, fps, As, fy, Asc, fy, Ff, fc, bw);
    fprintf('        = %.4f in\n\n', a);
    fprintf('    A_fl  = (b-bw)*hf  = (%.0f-%.2f)*%.1f = %.3f in^2  @ hf/2 = %.4f in from top\n', ...
            b, bw, hf, A_fl, hf/2);
    fprintf('    A_web = bw*a        = %.2f*%.4f    = %.3f in^2  @ a/2  = %.4f in from top\n', ...
            bw, a, A_web, a/2);
    fprintf('    ybar  = (A_fl*hf/2 + A_web*a/2) / (A_fl + A_web)\n');
    fprintf('          = (%.3f*%.4f + %.3f*%.4f) / (%.3f + %.3f)\n', ...
            A_fl, hf/2, A_web, a/2, A_fl, A_web);
    fprintf('          = %.4f in  from top fiber\n\n', ybar);
    fprintf('\n--- NOMINAL MOMENT Mn ---\n');
    fprintf('  Formula: Mn = Aps*fps*(dp-ybar) + As*fy*(d-ybar) - As''*fy*(d''-ybar)\n');
    fprintf('         [With As = %.4f in^2,  As'' = %.4f in^2]\n', As, Asc);
    fprintf('  Mn = %.4f*%.2f*(%.4f - %.4f)\n', Aps_total, fps, dp, ybar);
    if As > 0
        fprintf('     + %.4f*%.0f*(%.4f - %.4f)\n', As, fy, d_s, ybar);
    end
    if Asc > 0
        fprintf('     - %.4f*%.0f*(%.4f - %.4f)\n', Asc, fy, d_sc, ybar);
    end
    fprintf('     = %.1f kip-in  (%.1f kip-ft)\n', Mn, Mn/12);
end

%% ======================================================================
%%  STEP 7 — Ductility check (net tensile strain εt)
%%
%%  ACI 318-19 §9.3.3:  εt ≥ 0.004 (minimum for flexural members)
%%  φ = 0.90  when  εt ≥ 0.005  (tension-controlled)
%%  φ interpolated linearly between 0.65 and 0.90 when 0.004 ≤ εt < 0.005
%% ======================================================================
c     = a / beta1;
eps_t = (dp - c) / c * 0.003;   % net tensile strain at extreme tension steel

fprintf('\n--- DUCTILITY CHECK ---\n');
fprintf('  c   = a/β1 = %.4f/%.2f = %.4f in\n', a, beta1, c);
fprintf('  εt  = (dp−c)/c × 0.003\n');
fprintf('      = (%.4f−%.4f)/%.4f × 0.003\n', dp, c, c);
fprintf('      = %.5f\n', eps_t);

if eps_t >= 0.005
    phi = 0.90;
    ductility_status = 'Tension-controlled  (εt ≥ 0.005)  →  φ = 0.90 ✓';
elseif eps_t >= 0.004
    phi = 0.65 + (eps_t - 0.002) / (0.005 - 0.002) * (0.90 - 0.65);
    ductility_status = sprintf('Transition  (0.004 ≤ εt < 0.005)  →  φ = %.3f', phi);
else
    phi = NaN;
    ductility_status = 'FAIL — εt < 0.004  (over-reinforced, NOT permitted by ACI §9.3.3)';
    warning('eps_t = %.5f < 0.004: section is over-reinforced (ACI §9.3.3)', eps_t);
end
fprintf('  %s\n', ductility_status);

phi_Mn = phi * Mn;
fprintf('  φMn = %.2f × %.1f = %.1f kip-in  (%.1f kip-ft)\n', phi, Mn, phi_Mn, phi_Mn/12);

%% ======================================================================
%%  STEP 8 — Factored moment demand Mu along span
%% ======================================================================
x = beam.x;
L = beam.L;

% Load intensities
Ac = A_sec;
w_sw  = Ac * loads.concrete_density / 1000;   % kip/in (Ac*lb/in^3 / 1000 = kip/in)
w_SDL = 2 * 120 * (150/1728) / 1000;   % kip/in (2" topping, 150 pcf)
w_LL  = 0.42 / 12;                     % kip/in (0.42 kip/ft live load)

wu = 1.2 * (w_sw + w_SDL) + 1.6 * w_LL;   % factored UDL (kip/in)
Mu = wu .* x .* (L - x) / 2;              % factored moment diagram (kip-in)
Mu_max = max(Mu);

fprintf('\n--- FACTORED LOADS ---\n');
fprintf('  w_sw  = %.6f kip/in  (%.4f kip/ft)\n', w_sw,  w_sw*12);
fprintf('  w_SDL = %.6f kip/in  (%.4f kip/ft)\n', w_SDL, w_SDL*12);
fprintf('  w_LL  = %.6f kip/in  (%.4f kip/ft)\n', w_LL,  w_LL*12);
fprintf('  wu    = 1.2·(%.6f+%.6f) + 1.6·%.6f\n', w_sw, w_SDL, w_LL);
fprintf('        = %.6f kip/in  (%.4f kip/ft)\n', wu, wu*12);
fprintf('  Mu_max= wu·L²/8 = %.6f×%.0f²/8\n', wu, L);
fprintf('        = %.1f kip-in  (%.1f kip-ft)\n', Mu_max, Mu_max/12);

%% ======================================================================
%%  STEP 9 — Minimum flexural strength (anti-brittle failure)
%%
%%  ACI 318-19 §9.6.1.2:   φMn ≥ 1.2·Mcr
%%
%%  Mcr = (fr + fce) · Sb
%%    fr  = 7.5·√f'c_psi / 1000  (ksi)
%%    fce = Fe/Ac + Fe·e/Sb       (effective compressive stress at bottom fiber)
%% ======================================================================
fr  = materials.fr;              % ksi  (7.5√6000/1000 = 0.581 ksi)
fce = Fe/Ac + Fe*e_mid/Sb;       % ksi  (compression at bottom = positive)
Mcr = (fr + fce) * Sb;           % kip-in

fprintf('\n--- MINIMUM STRENGTH (Mcr) ---\n');
fprintf('  fr   = 7.5√(%.0f) / 1000 = %.4f ksi\n', fc*1000, fr);
fprintf('  Fe   = %.2f kips  (total effective prestress)\n', Fe);
fprintf('  fce  = Fe/Ac + Fe·e/Sb\n');
fprintf('       = %.2f/%.3f + %.2f×%.4f/%.1f\n', Fe, Ac, Fe, e_mid, Sb);
fprintf('       = %.4f ksi\n', fce);
fprintf('  Mcr  = (fr + fce)·Sb = (%.4f + %.4f)×%.1f\n', fr, fce, Sb);
fprintf('       = %.1f kip-in  (%.1f kip-ft)\n', Mcr, Mcr/12);
fprintf('  1.2·Mcr = %.1f kip-in  (%.1f kip-ft)\n', 1.2*Mcr, 1.2*Mcr/12);

%% ======================================================================
%%  STEP 10 — Design checks
%% ======================================================================
fprintf('\n=== DESIGN CHECKS ===\n');

% --- Check 1: φMn ≥ Mu_max ---
DC_strength = Mu_max / phi_Mn;
if phi_Mn >= Mu_max
    str_strength = 'PASS';
else
    str_strength = 'FAIL';
end
fprintf('\n  [1] Flexural Strength:  φMn ≥ Mu\n');
fprintf('      φMn    = %.1f kip-in  (%.1f kip-ft)\n', phi_Mn, phi_Mn/12);
fprintf('      Mu_max = %.1f kip-in  (%.1f kip-ft)\n', Mu_max, Mu_max/12);
fprintf('      D/C    = %.3f   →  %s\n', DC_strength, str_strength);

% --- Check 2: φMn ≥ 1.2·Mcr ---
DC_mcr = 1.2 * Mcr / phi_Mn;
if phi_Mn >= 1.2 * Mcr
    str_mcr = 'PASS';
else
    str_mcr = 'FAIL';
end
fprintf('\n  [2] Minimum Strength:   φMn ≥ 1.2·Mcr\n');
fprintf('      φMn    = %.1f kip-in\n', phi_Mn);
fprintf('      1.2Mcr = %.1f kip-in\n', 1.2*Mcr);
fprintf('      ratio  = %.3f   →  %s\n', DC_mcr, str_mcr);

% --- Ductility (already checked above) ---
if eps_t >= 0.004
    fprintf('\n  [3] Ductility:  εt = %.5f ≥ 0.004   →  PASS\n', eps_t);
else
    fprintf('\n  [3] Ductility:  εt = %.5f < 0.004   →  FAIL\n', eps_t);
end

%% ======================================================================
%%  STEP 11 — Plots
%% ======================================================================

% --- Figure 1: Mu vs φMn envelope ---
fig1 = figure('Name','Mu_vs_phiMn');
plot(x/12, Mu/12,      'r-',  'LineWidth', 2.0, 'DisplayName', 'M_u (factored)');
hold on;
yline(phi_Mn/12, 'b--', 'LineWidth', 1.8, 'DisplayName', sprintf('\\phiM_n = %.0f kip-ft', phi_Mn/12));
yline(1.2*Mcr/12, 'g-.', 'LineWidth', 1.5, 'DisplayName', sprintf('1.2M_{cr} = %.0f kip-ft', 1.2*Mcr/12));
hold off;
xlabel('x  (ft)');
ylabel('Moment  (kip-ft)');
title({'Factored Demand vs. Nominal Capacity', ...
       sprintf('\\phiM_n = %.0f kip-ft  |  M_{u,max} = %.0f kip-ft  |  D/C = %.2f', ...
               phi_Mn/12, Mu_max/12, DC_strength)});
legend('Location','north');
grid on;  box on;
set(gca,'FontSize',11);

% Annotate pass/fail
if phi_Mn >= Mu_max
    text(L/2/12, Mu_max/12 * 0.6, sprintf('\\color{blue}%s', str_strength), ...
         'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
else
    text(L/2/12, phi_Mn/12 * 1.05, sprintf('\\color{red}%s', str_strength), ...
         'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% --- Figure 2: Strain distribution at midspan ---
fig2 = figure('Name','Strain_at_Midspan');
eps_cu = 0.003;
eps_ce = 0;   % concrete strain at tendon level due to effective prestress
              % (approximate; for simplicity, show linear profile only)

% Strain at extreme compression fiber: eps_cu = 0.003 at top
% Linear profile: strain zero at NA (depth c from top), eps_t at dp
y_plot  = [0, dp];       % in (measured from extreme compression fiber = top)
eps_prf = [eps_cu, -eps_t];  % tension negative convention for clarity

y_fiber = linspace(0, h, 200);    % 0 at top, h at bottom
eps_fiber = interp1([0, c, dp], [eps_cu, 0, -eps_t], y_fiber, 'linear', 'extrap');

ax = gca;
fill([eps_fiber, 0, 0], [y_fiber, y_fiber(end), y_fiber(1)], ...
     [0.85 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on;
plot(eps_fiber, y_fiber, 'b-', 'LineWidth', 2);
xline(0, 'k--', 'LineWidth', 0.8);
plot(0, c, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');   % neutral axis
plot(-eps_t, dp, 'rs', 'MarkerSize', 9, 'MarkerFaceColor', 'r');  % tendon
yline(h - hf, 'g--', 'LineWidth', 1, 'Label', 'Flange soffit');
set(gca, 'YDir', 'reverse');   % top = y=0 (comp. fiber) at top of plot
xlabel('Strain  (in/in)');
ylabel('Distance from top fiber  (in)');
title({'Strain Distribution at Midspan (Ultimate)', ...
       sprintf('c = %.3f in  |  \\epsilon_t = %.5f  |  \\phi = %.2f', c, eps_t, phi)});
legend({'','Strain profile','Neutral axis','NA depth','Tendon (ε_t)'}, ...
       'Location','southeast');
grid on;  box on;
set(gca,'FontSize',11);

% --- Figure 3: Key dimension diagram ---
fig3 = figure('Name','Section_Capacity_Summary');
ax3 = axes;
hold on;
% Draw section outline (simplified)
patch(section.vertices(:,1), section.vertices(:,2), [0.9 0.9 0.9], ...
      'EdgeColor','k','LineWidth',1.5);
% Mark key depths (from bottom)
yline(yc, 'b--', 'LineWidth',1.5, 'Label', sprintf('yc = %.2f in', yc));
yline(y_tendon_mid, 'r-', 'LineWidth',2.0, 'Label', sprintf('y_{p} = %.2f in', y_tendon_mid));
% Neutral axis from bottom = h - c
yline(h - c, 'k:', 'LineWidth',1.5, 'Label', sprintf('NA  (c=%.3f in from top)', c));
% Compression block bottom from = h - a
yline(h - a, 'm--', 'LineWidth',1.2, 'Label', sprintf('a = %.3f in', a));
xlabel('Width  (in)');
ylabel('Height from bottom  (in)');
title({'Double-T Cross-Section — Ultimate Capacity', ...
       sprintf('f_p_s = %.1f ksi  |  a = %.4f in  |  φM_n = %.0f kip-ft', ...
               fps, a, phi_Mn/12)});
xlim([-70, 70]);  ylim([-1, 30]);
grid on;  box on;  axis equal;
set(gca,'FontSize',11);
hold off;

%% ======================================================================
%%  STEP 12 — Save figures
%% ======================================================================
fig_list = {fig1, fig2, fig3};
for k = 1:numel(fig_list)
    fig_k    = fig_list{k};
    fig_name = get(fig_k, 'Name');
    if isempty(fig_name);  fig_name = sprintf('Ultimate_Fig%d', k);  end
    saveas(fig_k, fullfile(output_dir, [fig_name '.png']));
    try
        savefig(fig_k, fullfile(output_dir, [fig_name '.fig']));
    catch
    end
    fprintf('  Saved: %s\n', fig_name);
end

%% ======================================================================
%%  STEP 13 — Write text report  (professor-grade format)
%% ======================================================================
rpt_file = fullfile(output_dir, 'Ultimate_Report_Tsection.txt');
fid = fopen(rpt_file, 'w');

fprintf(fid, '==================================================================\n');
fprintf(fid, '   CEE 530 PRESTRESSED CONCRETE  —  PROJECT 1  (PART 3)\n');
fprintf(fid, '   Double-T Beam: Ultimate Strength — Nominal Moment Capacity\n');
fprintf(fid, '==================================================================\n');
fprintf(fid, '   Student      : Chidchanok Pleesudjai (Fen)\n');
fprintf(fid, '   Date         : %s\n', datestr(now, 'yyyy-mm-dd'));
fprintf(fid, '   Active code  : %s  [toggle design_code in ultimateDesign.m]\n', design_code);
fprintf(fid, '==================================================================\n\n');


% -----------------------------------------------------------------------
%  PART 1 — GIVEN DATA
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 1 — GIVEN DATA\n');
fprintf(fid, '==================================================================\n\n');

fprintf(fid, '1.1  BEAM GEOMETRY\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  Span:  L = (20 + 24 + 20) ft x 12 in/ft = %.0f in  (%.1f ft)\n', L, L/12);
fprintf(fid, '\n');

fprintf(fid, '1.2  CROSS-SECTION  (Double-T / TT, symmetric)\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  Section origin : y = 0 at BOTTOM of stems, positive upward\n');
fprintf(fid, '  Top flange     : %.0f in wide x %.0f in thick   (y = %.0f to %.0f in)\n', b, hf, h-hf, h);
fprintf(fid, '  Two stems      : bw = %.2f in total (2 stems x 4.75 in avg)\n', bw);
fprintf(fid, '  Total depth    : h = %.1f in\n', h);
fprintf(fid, '\n');

fprintf(fid, '1.3  MATERIAL PROPERTIES\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  Concrete    : f''c  = %.1f ksi   (28-day strength)\n', fc);
fprintf(fid, '               f''ci = %.1f ksi   (strength at transfer)\n', materials.fci);
fprintf(fid, '  PS steel    : fpu  = %.0f ksi   (ultimate tensile strength)\n', fpu);
fprintf(fid, '               fpy  = %.1f ksi   (fpy/fpu = %.2f)\n', fpy, fpy/fpu);
fprintf(fid, '               Eps  = %.0f ksi\n', materials.Eps);
fprintf(fid, '  Mild steel  : fy   = %.0f ksi,  Es = %.0f ksi\n', materials.fy, materials.Es);
fprintf(fid, '               As   = %.4f in^2  (tension,     depth d  = %.4f in from top)\n', As,  d_s);
fprintf(fid, '               As'' = %.4f in^2  (compression, depth d'' = %.4f in from top)\n', Asc, d_sc);
fprintf(fid, '\n');

fprintf(fid, '1.4  PRESTRESS AND LOADS\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  Prestress losses : %.0f%%  -->  eta = 1 - 0.%.0f = %.2f\n', losses*100, losses*100, eta);
if ~isempty(n_strands_design)
    fprintf(fid, '  Design strands   : %d strands x %.3f in^2/strand\n', n_strands_design, Aps_per_strand);
    fprintf(fid, '    Aps = %d x %.3f = %.4f in^2  (total)\n', n_strands_design, Aps_per_strand, Aps_total);
    fprintf(fid, '    NOTE: inputData.m defines %d strands; %.3f in^2 overridden to %.0f strands for design.\n', ...
            n_tend, n_tend*Aps_per_strand, n_strands_design);
else
    fprintf(fid, '  Strands from inputData: %d x %.3f in^2 = %.4f in^2\n', n_tend, tendons{1}.Aps, Aps_total);
end
fprintf(fid, '  Tendon y at midspan : y_p = %.4f in  (from bottom, all strands at y=6 in at midspan)\n', y_tendon_mid);
fprintf(fid, '\n');
fprintf(fid, '  Self-weight (auto from section area):\n');
fprintf(fid, '    w_sw = Ac x rho_c\n');
fprintf(fid, '         = %.4f in^2 x (150/1728 lb/in^3) / 1000  [lb to kip]\n', Ac);
fprintf(fid, '         = %.6f kip/in  =  %.4f kip/ft\n', w_sw, w_sw*12);
fprintf(fid, '\n');
fprintf(fid, '  Superimposed dead load  (2-in concrete topping, 150 pcf):\n');
fprintf(fid, '    w_SDL = 2 in x %.0f in x (150/1728 lb/in^3) / 1000\n', b);
fprintf(fid, '          = %.6f kip/in  =  %.4f kip/ft\n', w_SDL, w_SDL*12);
fprintf(fid, '\n');
fprintf(fid, '  Live load:\n');
fprintf(fid, '    w_LL  = 0.42 kip/ft / 12 in/ft\n');
fprintf(fid, '          = %.6f kip/in  =  0.4200 kip/ft\n', w_LL);
fprintf(fid, '\n');

% -----------------------------------------------------------------------
%  PART 2 — SECTION PROPERTIES
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 2 — SECTION PROPERTIES  (Shoelace / Green''s Theorem)\n');
fprintf(fid, '==================================================================\n\n');

fprintf(fid, '2.1  CROSS-SECTIONAL AREA\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  Formula (Shoelace):\n');
fprintf(fid, '    A = (1/2) * |Sum_{i=1}^{n} (x_i * y_{i+1} - x_{i+1} * y_i)|\n\n');
fprintf(fid, '  A = Ac = %.4f in^2\n\n', Ac);

fprintf(fid, '2.2  CENTROID  (from bottom, y = 0 at stem soffit)\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  Formula:\n');
fprintf(fid, '    yc = [Sum (x_i*y_{i+1} - x_{i+1}*y_i)*(y_i + y_{i+1})] / (6*A)\n\n');
fprintf(fid, '  yc = %.4f in   from bottom\n', yc);
fprintf(fid, '  yt = h - yc  =  %.4f - %.4f  =  %.4f in   (centroid to TOP fiber)\n', h, yc, yt);
fprintf(fid, '  yb = yc      =  %.4f in             (centroid to BOTTOM fiber)\n\n', yb);

fprintf(fid, '2.3  MOMENT OF INERTIA  (centroidal x-axis)\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  Formula:\n');
fprintf(fid, '    Ic = (1/12)*|Sum (x_i*y_{i+1}-x_{i+1}*y_i)*(y_i^2+y_i*y_{i+1}+y_{i+1}^2)|\n\n');
fprintf(fid, '  Ic = %.2f in^4\n\n', Ic);

fprintf(fid, '2.4  SECTION MODULI\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  St = Ic / yt  =  %.2f / %.4f  =  %.4f in^3   (top)\n', Ic, yt, St);
fprintf(fid, '  Sb = Ic / yb  =  %.2f / %.4f  =  %.4f in^3   (bottom)\n\n', Ic, yb, Sb);

% -----------------------------------------------------------------------
%  PART 3 — ACI CODE PARAMETERS
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 3 — ACI CODE PARAMETERS  [%s]\n', design_code);
fprintf(fid, '==================================================================\n\n');

fprintf(fid, '3.1  STRESS BLOCK FACTOR  beta_1  (ACI 318-19 Table 22.2.2.4.3)\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  Formula:\n');
fprintf(fid, '    beta_1 = 0.85 - 0.05 * (f''c - 4)     for f''c > 4 ksi\n');
fprintf(fid, '    Clamped to [0.65, 0.85]\n\n');
fprintf(fid, '  beta_1 = 0.85 - 0.05 * (%.1f - 4)\n', fc);
fprintf(fid, '         = 0.85 - 0.05 * %.1f\n', fc-4);
fprintf(fid, '         = %.2f\n\n', beta1);

fprintf(fid, '3.2  PRESTRESS FACTOR  gamma_p  (ACI 318-19 Table 22.3.2.1)\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  fpy / fpu = %.1f / %.0f = %.3f\n', fpy, fpu, fpy/fpu);
if strcmp(design_code, 'CEE530')
    fprintf(fid, '  fpy/fpu >= 0.85  -->  gamma_p = 0.40  (stress-relieved strand) [CEE530]\n\n');
else
    fprintf(fid, '  fpy/fpu >= 0.90  -->  gamma_p = 0.28  (low-relaxation strand)  [ACI-318-19]\n\n');
end

fprintf(fid, '3.3  EFFECTIVE PRESTRESS CHECK  (ACI validity condition: fse >= 0.5*fpu)\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  fpi (initial) = Pi / Aps = 25.0 / %.3f = %.2f ksi\n', Aps_per_strand, 25.0/Aps_per_strand);
fprintf(fid, '  fse (effective after losses):\n');
fprintf(fid, '    fse = fpi * eta  =  %.2f * %.2f  =  %.2f ksi\n', fse_avg/eta, eta, fse_avg);
fprintf(fid, '    0.5 * fpu = 0.5 * %.0f = %.1f ksi\n', fpu, 0.5*fpu);
if fse_avg >= 0.5*fpu
    fprintf(fid, '    fse = %.2f ksi  >=  0.5*fpu = %.1f ksi   -->  VALID  ✓\n\n', fse_avg, 0.5*fpu);
else
    fprintf(fid, '    fse = %.2f ksi  <   0.5*fpu = %.1f ksi   -->  WARNING: ACI Eq. not valid\n\n', fse_avg, 0.5*fpu);
end

% -----------------------------------------------------------------------
%  PART 4 — TENDON EFFECTIVE DEPTH
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 4 — TENDON EFFECTIVE DEPTH  dp\n');
fprintf(fid, '==================================================================\n\n');

fprintf(fid, '  dp = distance from extreme compression fiber (top) to\n');
fprintf(fid, '       centroid of prestressing steel (tension side)\n\n');
fprintf(fid, '  Formula:\n');
fprintf(fid, '    dp = h - y_p\n');
fprintf(fid, '       where  y_p = force-weighted centroid of tendons from bottom\n\n');
fprintf(fid, '  At midspan (all 4 original tendons at y = 6 in):\n');
fprintf(fid, '    y_p = %.4f in  (from bottom)\n', y_tendon_mid);
fprintf(fid, '\n');
fprintf(fid, '  dp = h - y_p\n');
fprintf(fid, '     = %.4f - %.4f\n', h, y_tendon_mid);
fprintf(fid, '     = %.4f in\n\n', dp);
fprintf(fid, '  Equivalently:\n');
fprintf(fid, '    e_mid = yc - y_p  =  %.4f - %.4f  =  %.4f in  (eccentricity at midspan)\n', yc, y_tendon_mid, e_mid);
fprintf(fid, '    dp    = yt + e    =  %.4f + %.4f  =  %.4f in   ✓ (same result)\n\n', yt, e_mid, yt+e_mid);

% -----------------------------------------------------------------------
%  PART 5 — FACTORED MOMENT DEMAND  Mu
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 5 — FACTORED MOMENT DEMAND  Mu\n');
fprintf(fid, '==================================================================\n\n');

fprintf(fid, '  ACI 318-19 load combination (gravity only, Eq. 5.3.1b):\n');
fprintf(fid, '    U = 1.2D + 1.6L\n\n');
fprintf(fid, '  Factored uniform load:\n');
fprintf(fid, '    wu = 1.2 * (w_sw + w_SDL)  +  1.6 * w_LL\n');
fprintf(fid, '       = 1.2 * (%.6f + %.6f)  +  1.6 * %.6f\n', w_sw, w_SDL, w_LL);
fprintf(fid, '       = 1.2 * %.6f  +  %.6f\n', w_sw+w_SDL, 1.6*w_LL);
fprintf(fid, '       = %.6f  +  %.6f\n', 1.2*(w_sw+w_SDL), 1.6*w_LL);
fprintf(fid, '       = %.6f kip/in  =  %.4f kip/ft\n\n', wu, wu*12);
fprintf(fid, '  Maximum factored moment (at midspan, simply supported):\n');
fprintf(fid, '    Mu = wu * L^2 / 8\n');
fprintf(fid, '       = %.6f * %.0f^2 / 8\n', wu, L);
fprintf(fid, '       = %.6f * %.0f / 8\n', wu, L^2);
fprintf(fid, '       = %.2f kip-in  =  %.2f kip-ft\n\n', Mu_max, Mu_max/12);

% -----------------------------------------------------------------------
%  PART 6 — STRESS IN PRESTRESSING STEEL  fps
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 6 — STRESS IN PRESTRESSING STEEL  fps\n');
fprintf(fid, '==================================================================\n\n');

fprintf(fid, '  ACI 318-19 Eq. 22.3.2.1 (bonded tendons, general form):\n\n');
fprintf(fid, '    fps = fpu * [ 1 - (gamma_p / beta_1)\n');
fprintf(fid, '                    * ( rho_p * fpu/f''c\n');
fprintf(fid, '                      + (d / dp)  * omega\n');
fprintf(fid, '                      - (d''/ dp)  * omega'' ) ]\n\n');
fprintf(fid, '  Where:\n');
fprintf(fid, '    rho_p   = Aps / (b_eff * dp)         (prestress steel ratio)\n');
fprintf(fid, '    omega   = rho_s  * fy / f''c           (tension mild steel index)\n');
fprintf(fid, '    omega'' = rho_s'' * fy / f''c           (compression mild steel index)\n');
fprintf(fid, '    d       = depth to tension mild steel centroid from top (in)\n');
fprintf(fid, '    d''      = depth to compression mild steel centroid from top (in)\n\n');
fprintf(fid, '  Parameters:\n');
fprintf(fid, '    fpu     = %.0f ksi\n', fpu);
fprintf(fid, '    gamma_p = %.2f\n', gamma_p);
fprintf(fid, '    beta_1  = %.2f\n', beta1);
fprintf(fid, '    f''c     = %.1f ksi\n', fc);
fprintf(fid, '    fy      = %.0f ksi\n', fy);
fprintf(fid, '    As      = %.4f in^2   (tension mild steel)\n', As);
fprintf(fid, '    As''     = %.4f in^2   (compression mild steel)\n\n', Asc);
fprintf(fid, '  --- MILD STEEL SUBSTITUTION ---\n\n');
if As == 0 && Asc == 0
    fprintf(fid, '    As = 0   -->   omega  = rho_s  * fy/f''c = 0\n');
    fprintf(fid, '    As''= 0   -->   omega'' = rho_s'' * fy/f''c = 0\n\n');
    fprintf(fid, '    Substituting omega = 0,  omega'' = 0:\n\n');
    fprintf(fid, '    fps = fpu * [ 1 - (gamma_p / beta_1) * ( rho_p * fpu/f''c + 0 - 0 ) ]\n\n');
    fprintf(fid, '    Simplified to (prestress-only form):\n\n');
    fprintf(fid, '    fps = fpu * [ 1 - (gamma_p / beta_1) * (rho_p * fpu / f''c) ]\n\n');
else
    fprintf(fid, '    omega  = (As/(b_eff*dp))  * fy/f''c = (%.4f/(%.0f*%.4f)) * %.0f/%.1f = %.6f\n', ...
            As, b, dp, fy, fc, omega);
    fprintf(fid, '    omega'' = (As''/(b_eff*dp)) * fy/f''c = (%.4f/(%.0f*%.4f)) * %.0f/%.1f = %.6f\n\n', ...
            Asc, b, dp, fy, fc, omega_c);
    fprintf(fid, '    fps = fpu * [ 1 - (gamma_p/beta_1) * (rho_p*fpu/f''c + (d/dp)*omega - (d''/dp)*omega'') ]\n\n');
end

fprintf(fid, '  --- STEP A: Trial compression block (assume rectangular) ---\n\n');
fprintf(fid, '    Use full flange width:  b = %.0f in\n\n', b);
fprintf(fid, '    rho_p (trial) = Aps / (b * dp)\n');
fprintf(fid, '                  = %.4f / (%.0f * %.4f)\n', Aps_total, b, dp);
fprintf(fid, '                  = %.6f\n\n', rho_p_trial);
fprintf(fid, '    fps (trial) = %.0f * [ 1 - (%.2f / %.2f) * (%.6f * %.0f / %.1f) ]\n', ...
        fpu, gamma_p, beta1, rho_p_trial, fpu, fc);
fprintf(fid, '               = %.0f * [ 1 - %.4f * %.4f ]\n', ...
        fpu, gamma_p/beta1, rho_p_trial*fpu/fc);
fprintf(fid, '               = %.0f * [ 1 - %.6f ]\n', fpu, (gamma_p/beta1)*rho_p_trial*fpu/fc);
fprintf(fid, '               = %.2f ksi\n\n', fps_trial);
fprintf(fid, '    a_trial = Aps * fps_trial / (0.85 * f''c * b)\n');
fprintf(fid, '            = %.4f * %.2f / (0.85 * %.1f * %.0f)\n', Aps_total, fps_trial, fc, b);
fprintf(fid, '            = %.4f / %.4f\n', Aps_total*fps_trial, 0.85*fc*b);
fprintf(fid, '            = %.4f in\n\n', a_trial);

if strcmp(section_type(1), 'R')
    fprintf(fid, '    a_trial = %.4f in  <=  hf = %.1f in   -->  RECTANGULAR section  ✓\n', a_trial, hf);
    fprintf(fid, '    Compression block lies entirely within the flange.\n');
    fprintf(fid, '    Use b = %.0f in (full flange width) for all calculations.\n\n', b);
    fprintf(fid, '  --- STEP B: Final fps and compression depth  (rectangular) ---\n\n');
    fprintf(fid, '    rho_p = %.6f  (same as trial, b = %.0f in)\n\n', rho_p_trial, b);
    fprintf(fid, '    fps   = %.0f * [ 1 - (%.2f / %.2f) * (%.6f * %.0f / %.1f) ]\n', ...
            fpu, gamma_p, beta1, rho_p_trial, fpu, fc);
    fprintf(fid, '          = %.2f ksi\n\n', fps);
    fprintf(fid, '    a     = Aps * fps / (0.85 * f''c * b)\n');
    fprintf(fid, '          = %.4f * %.2f / (0.85 * %.1f * %.0f)\n', Aps_total, fps, fc, b);
    fprintf(fid, '          = %.4f in\n\n', a);
    fprintf(fid, '    Verify:  a = %.4f in  <=  hf = %.1f in   ✓\n\n', a, hf);
else
    fprintf(fid, '    a_trial = %.4f in  >  hf = %.1f in   -->  T-SECTION required\n', a_trial, hf);
    fprintf(fid, '    Neutral axis extends into the web. Must use bw and T-section formula.\n\n');
    fprintf(fid, '  --- STEP B: T-section recomputation ---\n\n');
    fprintf(fid, '    Use web width only:  bw = %.2f in\n\n', bw);
    fprintf(fid, '    rho_pw = Aps / (bw * dp)\n');
    fprintf(fid, '           = %.4f / (%.2f * %.4f)\n', Aps_total, bw, dp);
    fprintf(fid, '           = %.6f\n\n', rho_pw);
    fprintf(fid, '    fps    = %.0f * [ 1 - (%.2f / %.2f) * (%.6f * %.0f / %.1f) ]\n', ...
            fpu, gamma_p, beta1, rho_pw, fpu, fc);
    fprintf(fid, '           = %.2f ksi\n\n', fps);
    fprintf(fid, '    Compression depth (web portion only):\n');
    fprintf(fid, '    a = [Aps*fps - 0.85*f''c*(b-bw)*hf] / (0.85*f''c*bw)\n');
    fprintf(fid, '      = [%.4f*%.2f - 0.85*%.1f*(%.0f-%.2f)*%.1f] / (0.85*%.1f*%.2f)\n', ...
            Aps_total, fps, fc, b, bw, hf, fc, bw);
    fprintf(fid, '      = %.4f in\n\n', a);
    fprintf(fid, '    Centroid of T-shaped compression block from top fiber:\n');
    fprintf(fid, '      A_fl  = (b-bw) * hf         = (%.0f-%.2f) * %.1f = %.3f in^2\n', b, bw, hf, A_fl);
    fprintf(fid, '      A_web = bw * a               = %.2f * %.4f = %.3f in^2\n', bw, a, A_web);
    fprintf(fid, '      ybar  = (A_fl * hf/2  +  A_web * a/2) / (A_fl + A_web)\n');
    fprintf(fid, '            = (%.3f*%.4f + %.3f*%.4f) / (%.3f + %.3f)\n', ...
            A_fl, hf/2, A_web, a/2, A_fl, A_web);
    fprintf(fid, '            = %.4f in  from top fiber\n\n', ybar);
end

% -----------------------------------------------------------------------
%  PART 7 — NOMINAL MOMENT CAPACITY  Mn
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 7 — NOMINAL MOMENT CAPACITY  Mn\n');
fprintf(fid, '==================================================================\n\n');

if strcmp(section_type(1), 'R')
    fprintf(fid, '  RECTANGULAR SECTION  (a = %.4f in  <=  hf = %.1f in)\n\n', a, hf);
    fprintf(fid, '  Formula:\n');
    fprintf(fid, '    Mn = Aps * fps * (dp - a/2)\n\n');
    fprintf(fid, '  Substitution:\n');
    fprintf(fid, '    Mn = %.4f * %.2f * (%.4f - %.4f/2)\n', Aps_total, fps, dp, a);
    fprintf(fid, '       = %.4f * %.2f * (%.4f - %.4f)\n', Aps_total, fps, dp, a/2);
    fprintf(fid, '       = %.4f * %.4f\n', Aps_total*fps, dp-a/2);
    fprintf(fid, '       = %.2f kip-in\n', Mn);
    fprintf(fid, '       = %.2f kip-ft\n\n', Mn/12);
else
    fprintf(fid, '  T-SECTION  (a = %.4f in  >  hf = %.1f in)\n\n', a, hf);
    fprintf(fid, '  Formula (moment about centroid of T-block, not about a/2):\n');
    fprintf(fid, '    Mn = Aps * fps * (dp - ybar)\n\n');
    fprintf(fid, '  Substitution:\n');
    fprintf(fid, '    Mn = %.4f * %.2f * (%.4f - %.4f)\n', Aps_total, fps, dp, ybar);
    fprintf(fid, '       = %.4f * %.4f\n', Aps_total*fps, dp-ybar);
    fprintf(fid, '       = %.2f kip-in\n', Mn);
    fprintf(fid, '       = %.2f kip-ft\n\n', Mn/12);
end

% -----------------------------------------------------------------------
%  PART 8 — DUCTILITY CHECK  (Net Tensile Strain)
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 8 — DUCTILITY CHECK  (Net Tensile Strain  et)\n');
fprintf(fid, '==================================================================\n\n');

fprintf(fid, '  ACI 318-19 §9.3.3:\n');
fprintf(fid, '    Minimum  et >= 0.004  for flexural members\n');
fprintf(fid, '    et >= 0.005  -->  tension-controlled  -->  phi = 0.90\n');
fprintf(fid, '    0.004 <= et < 0.005  -->  transition zone  -->  phi interpolated\n\n');

fprintf(fid, '  Depth to neutral axis:\n');
fprintf(fid, '    Formula:  c = a / beta_1\n');
fprintf(fid, '    c = %.4f / %.2f\n', a, beta1);
fprintf(fid, '      = %.4f in\n\n', c);

fprintf(fid, '  Net tensile strain at extreme prestressing steel level:\n');
fprintf(fid, '    Formula:  et = (dp - c) / c * ecu      (ecu = 0.003)\n');
fprintf(fid, '    et = (%.4f - %.4f) / %.4f * 0.003\n', dp, c, c);
fprintf(fid, '       = %.4f / %.4f * 0.003\n', dp-c, c);
fprintf(fid, '       = %.5f\n\n', eps_t);

if eps_t >= 0.005
    fprintf(fid, '  et = %.5f  >=  0.005   -->  TENSION-CONTROLLED\n', eps_t);
    fprintf(fid, '  phi = 0.90   ✓\n\n');
elseif eps_t >= 0.004
    fprintf(fid, '  et = %.5f  (0.004 <= et < 0.005)  -->  TRANSITION ZONE\n', eps_t);
    fprintf(fid, '  phi = 0.65 + (et - 0.002) / (0.005 - 0.002) * (0.90 - 0.65)\n');
    fprintf(fid, '      = 0.65 + (%.5f - 0.002) / 0.003 * 0.25\n', eps_t);
    fprintf(fid, '      = %.3f\n\n', phi);
else
    fprintf(fid, '  et = %.5f  <  0.004   -->  *** OVER-REINFORCED — NOT PERMITTED ***\n\n', eps_t);
end

fprintf(fid, '  Reduced nominal moment capacity:\n');
fprintf(fid, '    phi*Mn = %.2f * %.2f\n', phi, Mn);
fprintf(fid, '           = %.2f kip-in  =  %.2f kip-ft\n\n', phi_Mn, phi_Mn/12);

% -----------------------------------------------------------------------
%  PART 9 — CRACKING MOMENT  Mcr
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 9 — CRACKING MOMENT  Mcr  (ACI 318-19 §9.6.1.2)\n');
fprintf(fid, '==================================================================\n\n');

fprintf(fid, '  Purpose: Ensure phi*Mn >= 1.2*Mcr (prevent brittle fracture at first crack)\n\n');
fprintf(fid, '  Formula:\n');
fprintf(fid, '    Mcr = (fr + fce) * Sb\n\n');
fprintf(fid, '  Modulus of rupture:\n');
fprintf(fid, '    fr = 7.5 * lambda * sqrt(f''c_psi) / 1000      (lambda = 1.0 normal weight)\n');
fprintf(fid, '       = 7.5 * 1.0 * sqrt(%.0f) / 1000\n', fc*1000);
fprintf(fid, '       = 7.5 * %.4f / 1000\n', sqrt(fc*1000));
fprintf(fid, '       = %.4f ksi  (tension, positive in this formula)\n\n', fr);
fprintf(fid, '  Effective prestress at BOTTOM fiber after losses:\n');
fprintf(fid, '    fce = Fe/Ac  +  Fe*e/Sb\n');
fprintf(fid, '    where:\n');
fprintf(fid, '      Fe  = Aps * fse  =  %.4f * %.2f  =  %.2f kips\n', Aps_total, fse_avg, Fe);
fprintf(fid, '      e   = e_mid      =  %.4f in  (eccentricity at midspan)\n', e_mid);
fprintf(fid, '      Ac  = %.4f in^2\n', Ac);
fprintf(fid, '      Sb  = %.4f in^3\n\n', Sb);
fprintf(fid, '    fce = %.2f / %.4f  +  %.2f * %.4f / %.4f\n', Fe, Ac, Fe, e_mid, Sb);
fprintf(fid, '        = %.4f  +  %.4f\n', Fe/Ac, Fe*e_mid/Sb);
fprintf(fid, '        = %.4f ksi  (compression at bottom fiber, positive)\n\n', fce);
fprintf(fid, '  Cracking moment:\n');
fprintf(fid, '    Mcr = (fr + fce) * Sb\n');
fprintf(fid, '        = (%.4f + %.4f) * %.4f\n', fr, fce, Sb);
fprintf(fid, '        = %.4f * %.4f\n', fr+fce, Sb);
fprintf(fid, '        = %.2f kip-in  =  %.2f kip-ft\n\n', Mcr, Mcr/12);
fprintf(fid, '  Minimum strength requirement:\n');
fprintf(fid, '    1.2 * Mcr = 1.2 * %.2f  =  %.2f kip-in  =  %.2f kip-ft\n\n', Mcr, 1.2*Mcr, 1.2*Mcr/12);

% -----------------------------------------------------------------------
%  PART 10 — DESIGN CHECKS SUMMARY
% -----------------------------------------------------------------------
fprintf(fid, '==================================================================\n');
fprintf(fid, '   PART 10 — DESIGN CHECKS SUMMARY\n');
fprintf(fid, '==================================================================\n\n');

% ---- Check 1: phi*Mn >= Mu ----
fprintf(fid, '10.1  FLEXURAL STRENGTH CHECK:  phi*Mn  >=  Mu\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  phi*Mn  = %.2f kip-in  (%.2f kip-ft)\n', phi_Mn, phi_Mn/12);
fprintf(fid, '  Mu_max  = %.2f kip-in  (%.2f kip-ft)\n', Mu_max, Mu_max/12);
fprintf(fid, '  D/C     = Mu / phi*Mn  =  %.2f / %.2f  =  %.3f\n', Mu_max, phi_Mn, DC_strength);
if strcmp(str_strength, 'PASS')
    fprintf(fid, '  RESULT:  phi*Mn = %.2f  >=  Mu = %.2f kip-in   -->  PASS  ✓\n\n', phi_Mn, Mu_max);
else
    fprintf(fid, '  RESULT:  phi*Mn = %.2f  <   Mu = %.2f kip-in   -->  FAIL  ✗\n\n', phi_Mn, Mu_max);
end

% ---- Check 2: phi*Mn >= 1.2*Mcr ----
fprintf(fid, '10.2  MINIMUM REINFORCEMENT CHECK:  phi*Mn  >=  1.2*Mcr\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  phi*Mn  = %.2f kip-in\n', phi_Mn);
fprintf(fid, '  1.2*Mcr = %.2f kip-in\n', 1.2*Mcr);
fprintf(fid, '  Ratio   = phi*Mn / (1.2*Mcr)  =  %.2f / %.2f  =  %.3f\n', phi_Mn, 1.2*Mcr, phi_Mn/(1.2*Mcr));
if strcmp(str_mcr, 'PASS')
    fprintf(fid, '  RESULT:  phi*Mn = %.2f  >=  1.2*Mcr = %.2f kip-in   -->  PASS  ✓\n\n', phi_Mn, 1.2*Mcr);
else
    fprintf(fid, '  RESULT:  phi*Mn = %.2f  <   1.2*Mcr = %.2f kip-in   -->  FAIL  ✗\n\n', phi_Mn, 1.2*Mcr);
end

% ---- Check 3: Ductility ----
if eps_t >= 0.004;  duct_str = 'PASS';  else;  duct_str = 'FAIL';  end
fprintf(fid, '10.3  DUCTILITY REQUIREMENT:  et  >=  0.004\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  et = %.5f\n', eps_t);
if strcmp(duct_str, 'PASS')
    fprintf(fid, '  RESULT:  et = %.5f  >=  0.004   -->  PASS  ✓\n\n', eps_t);
else
    fprintf(fid, '  RESULT:  et = %.5f  <   0.004   -->  FAIL  ✗  (over-reinforced)\n\n', eps_t);
end

% ---- Final summary table ----
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  FINAL SUMMARY TABLE\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  %-35s  %-20s  %-20s  %s\n', 'Check', 'Demand', 'Capacity', 'Status');
fprintf(fid, '  %-35s  %-20s  %-20s  %s\n', repmat('-',1,35), repmat('-',1,20), repmat('-',1,20), repmat('-',1,6));
fprintf(fid, '  %-35s  %-20s  %-20s  %s\n', ...
    'Flexural strength: phi*Mn >= Mu', ...
    sprintf('Mu = %.1f kip-in', Mu_max), ...
    sprintf('phi*Mn = %.1f kip-in', phi_Mn), str_strength);
fprintf(fid, '  %-35s  %-20s  %-20s  %s\n', ...
    'Min. reinforcement: phi*Mn >= 1.2Mcr', ...
    sprintf('1.2Mcr = %.1f kip-in', 1.2*Mcr), ...
    sprintf('phi*Mn = %.1f kip-in', phi_Mn), str_mcr);
fprintf(fid, '  %-35s  %-20s  %-20s  %s\n', ...
    'Ductility: et >= 0.004', ...
    sprintf('et = %.5f', eps_t), '>= 0.004', duct_str);
fprintf(fid, '\n');

% ---- Key results ----
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  KEY RESULTS\n');
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, '  Design code   : %s\n', code_label);
fprintf(fid, '  Section type  : %s\n', section_type);
if ~isempty(n_strands_design)
    n_used = n_strands_design;
else
    n_used = n_tend;
end
fprintf(fid, '  Aps           = %.4f in^2  (%d strands x %.3f in^2)\n', Aps_total, n_used, Aps_per_strand);
fprintf(fid, '  dp            = %.4f in\n', dp);
fprintf(fid, '  rho_p         = %.6f\n', rho_p_trial);
fprintf(fid, '  fps           = %.2f ksi\n', fps);
fprintf(fid, '  a             = %.4f in\n', a);
fprintf(fid, '  c             = %.4f in\n', c);
fprintf(fid, '  et            = %.5f\n', eps_t);
fprintf(fid, '  phi           = %.2f\n', phi);
fprintf(fid, '  Mn            = %.2f kip-in  (%.2f kip-ft)\n', Mn, Mn/12);
fprintf(fid, '  phi*Mn        = %.2f kip-in  (%.2f kip-ft)\n', phi_Mn, phi_Mn/12);
fprintf(fid, '  Mu_max        = %.2f kip-in  (%.2f kip-ft)\n', Mu_max, Mu_max/12);
fprintf(fid, '  D/C           = %.3f\n', DC_strength);
fprintf(fid, '  fr            = %.4f ksi\n', fr);
fprintf(fid, '  fce           = %.4f ksi\n', fce);
fprintf(fid, '  Mcr           = %.2f kip-in  (%.2f kip-ft)\n', Mcr, Mcr/12);
fprintf(fid, '  1.2*Mcr       = %.2f kip-in\n', 1.2*Mcr);
fprintf(fid, '\n');
fprintf(fid, '==================================================================\n');
fprintf(fid, '   Generated by ultimateDesign.m\n');
fprintf(fid, '   Active: %s\n', code_label);
fprintf(fid, '==================================================================\n');

fclose(fid);
fprintf('\n  Report written to: %s\n', rpt_file);
fprintf('\n========================================\n');
fprintf('  ULTIMATE DESIGN COMPLETE\n');
fprintf('========================================\n');
