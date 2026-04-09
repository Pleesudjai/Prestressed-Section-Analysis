%% SHEAR_DESIGN.m
% ACI 318-19 Shear Design for Prestressed Double-T Beam (Project 2)
% Method: Naaman Sections 5.7 — Empirical Formulas (Vci / Vcw)
%
% HOW TO RUN:
%   Set MATLAB current folder to 08_matlab_program_analysis/, then >> shearDesign
%
% EQUATIONS (from class notes + Naaman pp. 221-260):
%   Vci = 0.6*lambda*sqrt(f'c)*bw*dp + Vd + (Vi/Mmax)*Mcr
%         with Vci >= 5*lambda*sqrt(f'c)*bw*dp   (CEE530; ACI 318-19 uses 1.7)
%   Mcr = (I/yb)*(6*lambda*sqrt(f'c) + fce - fd)
%   Vcw = (3.5*lambda*sqrt(f'c) + 0.3*fpc)*bw*dp + Vp
%   Vc  = min(Vci, Vcw)
%   Vu <= phi*(Vc + Vs)  =>  Vs >= Vu/phi - Vc
%   Av/s = Vs / (fy * d)     [vertical stirrups, beta=45 deg]

clear; clc; close all;

fprintf('======================================================\n');
fprintf('   SHEAR DESIGN  —  ACI 318-19  /  Naaman Sec. 5.7  \n');
fprintf('   Project 2 : Prestressed Double-T Beam             \n');
fprintf('======================================================\n');

%% ================================================================
%%  1.  LOAD INPUT DATA
%% ================================================================
[beam, section, materials, prestress, reinforcement, loads] = inputData();

%% ================================================================
%%  2.  PARAMETERS
%% ================================================================
fc      = materials.fc;        % f'c  (ksi)
fci     = materials.fci;       % f'ci (ksi)
lambda  = 1.0;                 % normal-weight concrete
phi_v   = 0.75;                % ACI 318-19 phi for shear

% Stirrup — assume #3 U-stirrups (closed)
Av_bar  = 2 * 0.11;           % in²  (2 legs × 0.11 in²/leg)
fy_s    = 60.0;                % ksi  stirrup steel yield

fpu     = materials.fpu;       % ksi
L       = beam.L;              % in
x       = beam.x;              % 1×n array of positions along beam
n       = length(x);

sqrtfc  = sqrt(fc  * 1000) / 1000;   % sqrt(f'c_psi) expressed in ksi
sqrtfci = sqrt(fci * 1000) / 1000;

%% ================================================================
%%  3.  SECTION GEOMETRY
%% ================================================================
h   = max(section.vertices(:,2)) - min(section.vertices(:,2));  % total depth (in)
Ac  = section.A;
Ix  = section.Ix;
yc  = section.yc;   % centroid from BOTTOM (in)
yb  = section.yb;   % centroid to bottom fiber (in)
yt  = section.yt;   % centroid to top fiber (in)

% Web width bw — computed from section.vertices at minimum level (stem bottom)
% For Double-T: vertices at y_min give the narrowest cross-section (stem tips)
% Pair consecutive x-values to find each stem's width, then sum.
y_min_v = min(section.vertices(:,2));
x_at_ymin = sort(section.vertices(abs(section.vertices(:,2) - y_min_v) < 1e-6, 1));
bw = 0;
for kv = 1:2:length(x_at_ymin)-1
    bw = bw + (x_at_ymin(kv+1) - x_at_ymin(kv));
end

fprintf('\n------ SECTION PROPERTIES ------\n');
fprintf('  h   = %.2f in\n',  h);
fprintf('  Ac  = %.3f in^2\n',Ac);
fprintf('  Ix  = %.1f in^4\n',Ix);
fprintf('  yc  = %.4f in   (from bottom)\n', yc);
fprintf('  yb  = %.4f in\n', yb);
fprintf('  yt  = %.4f in\n', yt);
fprintf('  bw  = %.2f in   (min web width at stem bottom, from vertices)\n', bw);

%% ================================================================
%%  4.  TENDON PROPERTIES
%% ================================================================
n_tend    = length(prestress.tendons);
Aps_total = sum(cellfun(@(t) t.Aps, prestress.tendons));

% Resultant PS centroid y(x) from bottom
y_ps = zeros(1,n);
for i = 1:n_tend
    t = prestress.tendons{i};
    y_ps = y_ps + t.Aps * t.y;
end
y_ps = y_ps / Aps_total;    % in (from bottom)

% Effective depth dp = h - y_ps  (from compression face)
dp_raw = h - y_ps;
dp     = max(dp_raw, 0.8*h);   % ACI: dp >= 0.8h

% Eccentricity of resultant (positive = below centroid)
e_ps = yc - y_ps;

% Effective prestress force (constant along span for bonded tendons)
eta      = 1 - prestress.losses;
Pi_total = sum(cellfun(@(t) t.Aps * t.fpi, prestress.tendons));
Pe_val   = Pi_total * eta;     % kips  (scalar)
Pe       = Pe_val * ones(1,n); % kips  (vector)

fprintf('\n------ PRESTRESS ------\n');
fprintf('  Aps_total = %.4f in^2\n', Aps_total);
fprintf('  Pi_total  = %.2f kips\n', Pi_total);
fprintf('  Pe_total  = %.2f kips  (eta = %.2f)\n', Pe_val, eta);
fprintf('  dp at midspan = %.2f in  (0.8h = %.2f in)\n', ...
        dp(round(n/2)), 0.8*h);
fprintf('  e at midspan  = %.3f in  (+ = below centroid)\n', e_ps(round(n/2)));

%% ================================================================
%%  5.  VERTICAL COMPONENT OF PRESTRESS (Vp)
%% ================================================================
% Computed analytically from y_profile of each tendon
% Vp > 0 = upward component (helps resist applied shear)
Vp = zeros(1,n);

for i = 1:n_tend
    t = prestress.tendons{i};
    if isfield(t,'y_profile') && ~isempty(t.y_profile)
        % Piecewise-linear profile → slope at each x
        dy_dx = zeros(1,n);
        xp = t.x_profile;
        yp = t.y_profile;
        for seg = 1:length(xp)-1
            mask = (x >= xp(seg)) & (x <= xp(seg+1));
            dy_dx(mask) = (yp(seg+1)-yp(seg)) / (xp(seg+1)-xp(seg));
        end
        theta_i = atan(abs(dy_dx));
        Pe_i    = t.Aps * t.fpi * eta;
        % Left side: downward slope of tendon → upward Vp
        % Right side: upward slope → also upward Vp (symmetric)
        Vp = Vp + Pe_i * sin(theta_i);
    end
end

fprintf('  Vp at support = %.3f kips  (upward, harped tendons)\n', Vp(1));
fprintf('  Vp at midspan = %.3f kips\n', Vp(round(n/2)));

%% ================================================================
%%  6.  APPLIED LOADS
%% ================================================================
% Self-weight
w_sw  = Ac * (150/1728) / 1000;    % kip/in  [150 pcf = 150/1728 lb/in^3]

% SDL (from loads.distributed)
w_sdl = 0;
if isfield(loads,'distributed') && ~isempty(loads.distributed)
    for k = 1:size(loads.distributed,1)
        w_sdl = w_sdl + abs(loads.distributed(k,3));
    end
end

% Live load
w_ll = 0;
if isfield(loads,'distributed_live') && ~isempty(loads.distributed_live)
    for k = 1:size(loads.distributed_live,1)
        w_ll = w_ll + abs(loads.distributed_live(k,3));
    end
end

w_DL = w_sw + w_sdl;                   % total unfactored DL (kip/in)
w_u  = 1.2 * w_DL + 1.6 * w_ll;       % total factored load (kip/in)

fprintf('\n------ LOADS ------\n');
fprintf('  w_sw  = %.5f kip/in  = %.4f kip/ft\n', w_sw,  w_sw*12);
fprintf('  w_sdl = %.5f kip/in  = %.4f kip/ft\n', w_sdl, w_sdl*12);
fprintf('  w_ll  = %.5f kip/in  = %.4f kip/ft\n', w_ll,  w_ll*12);
fprintf('  w_DL  = %.5f kip/in  = %.4f kip/ft  [SW+SDL]\n', w_DL, w_DL*12);
fprintf('  w_u   = %.5f kip/in  = %.4f kip/ft  [1.2DL+1.6LL]\n', w_u, w_u*12);

%% ================================================================
%%  7.  INTERNAL FORCES  (simply supported beam)
%% ================================================================
% Unfactored DL shear and moment (for Vd, fd, Mcr)
Vd_signed = w_DL * (L/2 - x);          % kips (+ left, - right)
Md        = w_DL * x .* (L - x) / 2;  % kip-in

% Factored total shear |Vu| and moment Mu
Vu_abs = abs(w_u * (L/2 - x));         % kips
Mu     = w_u * x .* (L - x) / 2;      % kip-in

% Vi and Mmax for Vci:
%   Per notes: factored shear/moment from EXTERNAL loads only (SDL + LL), not SW
%   w_ext_factored = 1.2*w_sdl + 1.6*w_ll
w_ext = 1.2 * w_sdl + 1.6 * w_ll;
Vi    = abs(w_ext * (L/2 - x));         % kips
Mmax  = w_ext * x .* (L - x) / 2;     % kip-in
Mmax(Mmax < 1e-6) = 1e-6;              % guard against x=0,L singularity

fprintf('\n  Max Vu (at support) = %.2f kips\n', max(Vu_abs));
fprintf('  Max Mu (at midspan) = %.2f kip-in  = %.2f kip-ft\n', ...
        max(Mu), max(Mu)/12);

%% ================================================================
%%  8.  PRESTRESS STRESSES AT EXTREME BOTTOM FIBER
%% ================================================================
% fce = compressive stress at bottom fiber due to EFFECTIVE prestress (ksi, + = compression)
%       fce = Pe/Ac + Pe*e*yb/Ix
fce = Pe./Ac + Pe .* e_ps * yb ./ Ix;

% fd  = flexural stress at bottom fiber due to UNFACTORED DEAD LOAD
%       (positive = tensile at bottom, consistent with Naaman Mcr formula)
fd_bot = abs(Md) * yb ./ Ix;   % ksi

fprintf('\n------ STRESSES AT BOTTOM FIBER (midspan) ------\n');
mid = round(n/2);
fprintf('  fce  = %.4f ksi  (prestress compression)\n', fce(mid));
fprintf('  fd   = %.4f ksi  (DL flexure, tensile magnitude)\n', fd_bot(mid));

%% ================================================================
%%  9.  CRACKING MOMENT Mcr
%%      Mcr = (I/yb) * (6*lambda*sqrt(f'c) + fce - fd)
%% ================================================================
fr  = 6 * lambda * sqrtfc;    % modulus of rupture for shear  (ksi)
% Note: 6*sqrt(f'c_psi)/1000 = 6*sqrt(6000)/1000 = 0.4648 ksi

Mcr_raw = (Ix / yb) * (fr + fce - fd_bot);   % kip-in (may be negative)
Mcr = max(Mcr_raw, 0);

fprintf('\n------ CRACKING MOMENT (midspan) ------\n');
fprintf('  fr   = 6*lambda*sqrt(f''c) = %.4f ksi\n', fr);
fprintf('  Mcr  = (I/yb)*(fr + fce - fd) = %.1f kip-in  = %.1f kip-ft\n', ...
        Mcr(mid), Mcr(mid)/12);
if Mcr_raw(mid) < 0
    fprintf('  NOTE: Mcr < 0 at midspan — DL tensile stress (fd=%.3f ksi) exceeds\n', fd_bot(mid));
    fprintf('        prestress + rupture (fce+fr=%.3f ksi). Section under-prestressed.\n', fce(mid)+fr);
    fprintf('        Vci = Vci_min governs at midspan. Vcw provides adequate capacity.\n');
end

%% ================================================================
%%  10. Vci  —  FLEXURAL-SHEAR STRENGTH
%%      Vci = 0.6*lambda*sqrt(f'c)*bw*dp + Vd + (Vi/Mmax)*Mcr
%%
%%  BOUNDS (ACI 318-19 Table 22.5.8.2):
%%      Vci >= Vci_min = 1.7 * lambda * sqrt(f'c) * bw * dp   [hard floor]
%%      Vci <= Vci_max = 5.0 * lambda * sqrt(f'c) * bw * dp   [CEE 530 cap]
%%
%%      ACI 318-19 hard floor (1.7):  prevents Vci from dropping
%%        unrealistically low near midspan where Vi→0 and Mcr→0.
%%      CEE 530 upper cap (5.0, Prof. Mobasher):  caps Vci so it cannot
%%        exceed the simplified-method upper bound.
%%
%%  WHY THE FLOOR EXISTS:
%%    Near midspan, Vi→0 while Mmax is large, so the third term vanishes.
%%    If Mcr <= 0 (DL tension > prestress + rupture), term 3 = 0 as well.
%%    Without a floor, Vci would drop to just Vd + 0.6*sqrt*bw*dp — an
%%    unrealistically low value.  The 1.7 floor recognizes that even
%%    after flexural cracking, concrete retains shear capacity from
%%    aggregate interlock, dowel action, and the compression zone.
%%
%%  WHY THE CAP EXISTS (CEE 530):
%%    Near supports, Vi/Mmax → ∞ which inflates Vci unrealistically.
%%    The 5.0 cap (equal to the simplified-method upper bound from
%%    ACI 318-19 Sec. 22.5.8.1) prevents Vci from exceeding what
%%    the simplified method would allow.  This is a course convention.
%% ================================================================
Vci_t1  = 0.6 * lambda * sqrtfc * bw * dp;   % concrete base
Vci_t2  = abs(Vd_signed);                    % unfactored DL shear
Vci_t3  = Vi .* Mcr ./ Mmax;                 % moment-shear interaction
% At supports Mmax→0: Vci→∞ → Vcw governs. Cap to avoid Inf in display.
Vci_t3(~isfinite(Vci_t3)) = 1e9;

% Full three-term Vci
Vci     = Vci_t1 + Vci_t2 + Vci_t3;

% Bounds — floor and cap:
%   Floor (ACI 318-19): Vci >= 1.7 * lambda * sqrt(f'c) * bw * dp
%   Cap   (CEE 530):    Vci <= 5.0 * lambda * sqrt(f'c) * bw * dp
Vci_min_ACI = 1.7 * lambda * sqrtfc * bw * dp;   % ACI 318-19 hard floor
Vci_max_CEE = 5.0 * lambda * sqrtfc * bw * dp;   % CEE 530 upper cap
Vci     = max(Vci, Vci_min_ACI);   % enforce floor
Vci     = min(Vci, Vci_max_CEE);   % enforce cap

fprintf('\n------ Vci at Critical Section & Midspan ------\n');
[~,ic] = min(abs(x - max(dp(1), 0.8*h)));    % index near critical section
fprintf('  [x = %.1f in = %.2f ft]:\n', x(ic), x(ic)/12);
fprintf('    term1 = 0.6*lambda*sqrt(f''c)*bw*dp = %.2f kips\n', Vci_t1(ic));
fprintf('    term2 = Vd                          = %.2f kips\n', Vci_t2(ic));
fprintf('    term3 = (Vi/Mmax)*Mcr               = %.2f kips\n', Vci_t3(ic));
fprintf('    Vci   =                               %.2f kips\n', Vci(ic));
fprintf('    Vci_min (ACI 318-19)  = 1.7*lambda*sqrt(f''c)*bw*dp  = %.2f kips  [floor]\n', Vci_min_ACI(ic));
fprintf('    Vci_max (CEE 530)     = 5.0*lambda*sqrt(f''c)*bw*dp  = %.2f kips  [cap]\n', Vci_max_CEE(ic));

%% ================================================================
%%  11. Vcw  —  WEB-SHEAR STRENGTH
%%      Vcw = (3.5*lambda*sqrt(f'c) + 0.3*fpc)*bw*dp + Vp
%%      fpc = Pe/Ac  (compressive stress at centroid; e-term = 0 at centroid)
%%      Note: if centroid is in FLANGE, use stress at flange/web interface.
%%            For this Double-T: yc is in the STEM (yc < 26"), use centroid.
%% ================================================================
fpc = Pe ./ Ac;    % ksi (compressive at centroid)

Vcw = (3.5 * lambda * sqrtfc + 0.3 * fpc) .* bw .* dp + Vp;

fprintf('\n------ Vcw at Critical Section ------\n');
fprintf('  [x = %.1f in]:\n', x(ic));
fprintf('    fpc = Pe/Ac       = %.4f ksi\n', fpc(ic));
fprintf('    (3.5*sqrt+0.3fpc)*bw*dp = %.2f kips\n', ...
        (3.5*lambda*sqrtfc + 0.3*fpc(ic)) * bw * dp(ic));
fprintf('    Vp                = %.3f kips\n', Vp(ic));
fprintf('    Vcw               = %.2f kips\n', Vcw(ic));

%% ================================================================
%%  12. Vc = min(Vci, Vcw)
%% ================================================================
Vc = min(Vci, Vcw);

fprintf('\n------ Vc = min(Vci, Vcw) at Critical Section ------\n');
fprintf('  Vci = %.2f kips | Vcw = %.2f kips → Vc = %.2f kips\n', ...
        Vci(ic), Vcw(ic), Vc(ic));
fprintf('  phi*Vc = %.2f kips\n', phi_v*Vc(ic));

%% ================================================================
%%  13. STIRRUP DESIGN
%%      Vs_req = Vu/phi - Vc  (>= 0)
%%      Av/s   = Vs / (fy * dp)   for vertical stirrups (beta = 45 deg)
%% ================================================================
Vs_req = max(Vu_abs / phi_v - Vc, 0);    % kips

% Required Av/s
Av_s_demand = Vs_req ./ (fy_s * dp);     % in^2/in

% Minimum Av/s  (ACI 318-19 Table 9.6.3.3 — prestressed member)
%   Criterion 1: 0.75*sqrt(f'c)*bw/fy
%   Criterion 2: 50*bw/fy  (in psi → ksi: 50/1000)
%   Criterion 3: Aps*fpu/(80*fy*d) * sqrt(d/bw)
Av_s_c1 = 0.75 * sqrtfc * bw / fy_s;           % ksi system
Av_s_c2 = 50/1000 * bw / fy_s;                  % 50 psi → ksi
Av_s_c3 = Aps_total * fpu ./ (80 * fy_s * dp) .* sqrt(dp / bw);
Av_s_min = max([Av_s_c1, Av_s_c2, max(Av_s_c3)]);

Av_s_govern = max(Av_s_demand, Av_s_min);       % governing Av/s

% Required spacing for selected stirrup bar
s_req = Av_bar ./ Av_s_govern;   % in

% Spacing limits
s_lim_basic  = min(3*h/4, 24);   % 3h/4, 24" max
s_lim_tight  = min(3*h/8, 12);   % 3h/8, 12" — when Vs > 4*sqrt(f'c)*bw*d
Vs_thresh    = 4 * lambda * sqrtfc * bw * dp;
s_max        = s_lim_basic * ones(1,n);
s_max(Vs_req > Vs_thresh) = s_lim_tight;

% Vs max check
Vs_max_allow = 8 * lambda * sqrtfc * bw * dp;

s_design = min(s_req, s_max);    % governed by demand OR spacing limit
s_design = max(floor(s_design),1); % round down, min 1"

% Flag overstress
if any(Vs_req > Vs_max_allow)
    fprintf('\n*** WARNING: Vs > 8*sqrt(f''c)*bw*d at some locations — ENLARGE SECTION ***\n');
else
    fprintf('\n  Vs_max check OK everywhere (Vs <= 8*sqrt(f''c)*bw*d = %.1f kips)\n', ...
            Vs_max_allow(ic));
end

%% ================================================================
%%  14. PRINT FULL CALCULATION SUMMARY
%% ================================================================
fprintf('\n');
fprintf('====================================================\n');
fprintf('   SHEAR DESIGN CALCULATION SUMMARY\n');
fprintf('====================================================\n');
fprintf('Material:  f''c = %.1f ksi  |  lambda = %.1f\n', fc, lambda);
fprintf('Section:   bw  = %.2f in   |  h = %.2f in  |  dp >= %.2f in\n', bw, h, 0.8*h);
fprintf('Loads:     wu  = %.5f kip/in  (1.2DL+1.6LL)\n', w_u);
fprintf('Stirrups:  #3 U-stirrup  Av = %.2f in^2  |  fy = %.0f ksi\n', Av_bar, fy_s);
fprintf('\n');
fprintf('  Minimum Av/s:\n');
fprintf('    Criterion 1 (0.75*sqrt(f''c)*bw/fy)            = %.5f in^2/in\n', Av_s_c1);
fprintf('    Criterion 2 (50*bw/fy)                         = %.5f in^2/in\n', Av_s_c2);
fprintf('    Criterion 3 (Aps*fpu/(80*fy*d)*sqrt(d/bw))    = %.5f in^2/in\n', max(Av_s_c3));
fprintf('    GOVERNING Av/s_min                             = %.5f in^2/in\n', Av_s_min);
fprintf('\n');
fprintf('  Spacing limits:\n');
fprintf('    Basic:  s_max = min(3h/4, 24") = %.1f"\n', s_lim_basic);
fprintf('    Tight:  s_max = min(3h/8, 12") = %.1f"  [when Vs > 4*sqrt(f''c)*bw*d]\n', s_lim_tight);
fprintf('\n');

fprintf('%-10s %-8s %-8s %-8s %-8s %-8s %-8s %-10s %-8s\n', ...
    'x (ft)','Vu(k)','Vci(k)','Vcw(k)','Vc(k)','phiVc','Vs_req','Av/s','s_des"');
fprintf('%s\n', repmat('-',1,82));

% Print at design sections from beam.x_plot_fractions
x_print = beam.x_plot_fractions * L;   % [0, 120, 240, 384] in
for k = 1:length(x_print)
    [~,ip] = min(abs(x - x_print(k)));
    % Display Vci — cap at 999 for readability (huge = Vcw governs)
    Vci_disp = min(Vci(ip), 999.9);
    if Vci(ip) > 999
        Vci_str = ' >999  ';
    else
        Vci_str = sprintf('%-8.2f', Vci_disp);
    end
    fprintf('%-10.1f %-8.2f %s %-8.2f %-8.2f %-8.2f %-8.2f %-10.5f %-8d\n', ...
        x(ip)/12, Vu_abs(ip), Vci_str, Vcw(ip), Vc(ip), phi_v*Vc(ip), ...
        Vs_req(ip), Av_s_govern(ip), round(s_design(ip)));
end

%% ================================================================
%%  15. DETAILED CHECK AT EACH DESIGN SECTION
%% ================================================================
x_design = beam.x_plot_fractions * L;   % [0, 120, 240, 384] in
for ks = 1:length(x_design)
    [~,ip] = min(abs(x - x_design(ks)));

    fprintf('\n====================================================\n');
    fprintf('   SECTION %d :  x = %.2f in = %.2f ft\n', ks, x(ip), x(ip)/12);
    fprintf('====================================================\n');
    fprintf('  dp          = %.3f in   (>= 0.8h = %.3f in)\n', dp(ip), 0.8*h);
    fprintf('  y_ps        = %.3f in   (from bottom)\n', y_ps(ip));
    fprintf('  e           = %.3f in   (+ = below centroid)\n', e_ps(ip));
    fprintf('  Vu          = %.3f kips\n', Vu_abs(ip));

    fprintf('\n  Mcr calculation:\n');
    fprintf('    fr  = 6*lambda*sqrt(%.0f psi) = %.4f ksi\n', fc*1000, fr);
    fprintf('    fce = Pe/Ac + Pe*e*yb/Ix     = %.4f ksi\n', fce(ip));
    fprintf('    fd  = Md*yb/Ix               = %.4f ksi\n', fd_bot(ip));
    fprintf('    fr+fce-fd                     = %.4f ksi\n', fr + fce(ip) - fd_bot(ip));
    if Mcr(ip) > 0
        fprintf('    Mcr = (I/yb)*(fr+fce-fd)     = %.2f kip-in\n', Mcr(ip));
    else
        fprintf('    Mcr = 0  (section cracked under DL alone)\n');
    end

    fprintf('\n  Vci:\n');
    fprintf('    0.6*lambda*sqrt(f''c)*bw*dp  = %.3f kips\n', Vci_t1(ip));
    fprintf('    Vd                          = %.3f kips\n', Vci_t2(ip));
    if Mcr(ip) > 0 && isfinite(Vci_t3(ip))
        fprintf('    (Vi/Mmax)*Mcr               = %.3f kips\n', Vci_t3(ip));
        fprintf('    Vci (computed)              = %.3f kips\n', Vci_t1(ip)+Vci_t2(ip)+Vci_t3(ip));
    else
        fprintf('    (Vi/Mmax)*Mcr               = 0 (Mcr=0 or Mmax->0)\n');
    end
    fprintf('    Vci_min (ACI, 1.7)           = %.3f kips  [floor]\n', Vci_min_ACI(ip));
    fprintf('    Vci_max (CEE 530, 5.0)       = %.3f kips  [cap]\n', Vci_max_CEE(ip));
    fprintf('    Vci (governing)             = %.3f kips\n', Vci(ip));

    fprintf('\n  Vcw:\n');
    fprintf('    fpc = Pe/Ac                 = %.4f ksi\n', fpc(ip));
    fprintf('    (3.5*sqrt+0.3fpc)*bw*dp     = %.3f kips\n', ...
            (3.5*lambda*sqrtfc+0.3*fpc(ip))*bw*dp(ip));
    fprintf('    Vp (harped tendons)         = %.3f kips\n', Vp(ip));
    fprintf('    Vcw                         = %.3f kips\n', Vcw(ip));

    fprintf('\n  Vc = min(Vci=%.2f, Vcw=%.2f) = %.2f kips\n', Vci(ip), Vcw(ip), Vc(ip));
    fprintf('  phi*Vc = %.2f kips\n', phi_v*Vc(ip));
    if Vu_abs(ip) <= phi_v * Vc(ip)
        fprintf('  Vu <= phi*Vc  -> Only minimum stirrups required at this section\n');
    else
        fprintf('  Vu > phi*Vc  -> Stirrups required\n');
    end
    fprintf('  Vs_req = max(Vu/phi - Vc, 0) = %.3f kips\n', Vs_req(ip));

    fprintf('\n  Stirrup design (#3 U-stirrups, Av = %.2f in^2):\n', Av_bar);
    fprintf('    Av/s_demand  = Vs/(fy*d)    = %.5f in^2/in\n', Av_s_demand(ip));
    fprintf('    Av/s_min     (governing)    = %.5f in^2/in\n', Av_s_min);
    fprintf('    Av/s_used    = max of above = %.5f in^2/in\n', Av_s_govern(ip));
    fprintf('    s_req  = Av / (Av/s)        = %.2f in\n', Av_bar/Av_s_govern(ip));
    fprintf('    s_max  = %.1f in\n', s_max(ip));
    fprintf('    s_design = min(s_req, s_max) = %d in   <- USE THIS\n', round(s_design(ip)));
end

%% ================================================================
%%  16. PLOTS
%% ================================================================
figure('Name','Shear Design Summary','Position',[30 30 1300 850]);

% Design section locations for vertical markers
x_des_ft = x_design / 12;

% Signed factored shear (positive at left, negative at right)
Vu_signed = w_u * (L/2 - x);   % kips (+ left half, - right half)

% Display cap: clip Vci for plotting (Vci→∞ near supports distorts scale)
y_cap = max(Vu_abs) * 1.5;   % cap at 1.5× max factored shear
phi_Vci_plot = min(phi_v * Vci, y_cap);
phi_Vcw_plot = phi_v * Vcw;

% ---- (a) Shear demand (signed) vs. capacity envelope ----
ax1 = subplot(2,2,1);
% Signed Vu — straight line from + to -
plot(x/12, Vu_signed, 'r-', 'LineWidth', 2.0, 'DisplayName', 'V_u');
hold on;
% Positive capacity envelope
plot(x/12, phi_v*Vc,       'b-',  'LineWidth', 2.0, 'DisplayName', '+\phi V_c');
plot(x/12, phi_Vci_plot,   'b--', 'LineWidth', 1.2, 'DisplayName', '+\phi V_{ci}');
plot(x/12, phi_Vcw_plot,   'g--', 'LineWidth', 1.2, 'DisplayName', '+\phi V_{cw}');
% Negative capacity envelope (mirror)
plot(x/12, -phi_v*Vc,      'b-',  'LineWidth', 2.0, 'HandleVisibility','off');
plot(x/12, -phi_Vci_plot,  'b--', 'LineWidth', 1.2, 'HandleVisibility','off');
plot(x/12, -phi_Vcw_plot,  'g--', 'LineWidth', 1.2, 'HandleVisibility','off');
% Zero line
plot([0 L/12], [0 0], 'k-', 'LineWidth', 0.5, 'HandleVisibility','off');
for ks = 1:length(x_des_ft)
    xline(x_des_ft(ks), 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
end
ylim([-y_cap, y_cap]);
xlabel('x  (ft)'); ylabel('Shear  (kips)');
title('(a)  Shear Demand vs. Concrete Capacity');
legend('Location','northeast','FontSize',8); grid on;

% ---- (b) Vs required ----
ax2 = subplot(2,2,2);
plot(x/12, Vs_req, 'r-', 'LineWidth', 2.0, 'DisplayName', 'V_s required');
hold on;
plot(x/12, Vs_max_allow, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', '8\lambda\surd f''_c b_w d  (max V_s)');
for ks = 1:length(x_des_ft)
    xline(x_des_ft(ks), 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
end
ylim([0, max(Vs_max_allow)*1.2]);
xlabel('x  (ft)'); ylabel('V_s  (kips)');
title('(b)  Required Stirrup Contribution V_s');
legend('Location','northeast','FontSize',8); grid on;

% ---- (c) Av/s ----
ax3 = subplot(2,2,3);
% Clip demand for display (spikes near supports where Mmax→0)
Av_s_demand_plot = min(Av_s_demand, max(Av_s_govern)*1.5);
plot(x/12, Av_s_govern, 'r-', 'LineWidth', 2.0, 'DisplayName', 'A_v/s  (governing)');
hold on;
plot(x/12, Av_s_demand_plot, 'm--', 'LineWidth', 1.2, 'DisplayName', 'A_v/s  (demand)');
yline(Av_s_min, 'k--', 'LineWidth', 1.5, 'DisplayName', 'A_v/s  (min)');
for ks = 1:length(x_des_ft)
    xline(x_des_ft(ks), 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
end
ylim([0, max(Av_s_govern)*1.5]);
xlabel('x  (ft)'); ylabel('A_v/s  (in^2/in)');
title('(c)  Required A_v/s   [#3 U-stirrups: A_v = 0.22 in^2]');
legend('Location','northeast','FontSize',8); grid on;

% ---- (d) Stirrup spacing ----
ax4 = subplot(2,2,4);
stairs(x/12, s_design, 'b-', 'LineWidth', 2.0, 'DisplayName', 's_{design}  (in)');
hold on;
plot(x/12, s_max, 'k--', 'LineWidth', 1.5, 'DisplayName', 's_{max}');
for ks = 1:length(x_des_ft)
    xline(x_des_ft(ks), 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
end
ylim([0, s_lim_basic + 4]);
xlabel('x  (ft)'); ylabel('Stirrup Spacing  s  (in)');
title('(d)  Stirrup Spacing  (#3 U-stirrups)');
legend('Location','southeast','FontSize',8); grid on;

sgtitle(sprintf('Shear Design — Double-T Beam   |   f''c = %.0f ksi,  bw = %.1f in,  h = %.0f in,  L = %.0f ft', ...
    fc, bw, h, L/12), 'FontSize', 13, 'FontWeight', 'bold');

% ---- Save ----
out_dir = fullfile('project_Project2', 'output');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
saveas(gcf, fullfile(out_dir,'ShearDesign.png'));
savefig(gcf, fullfile(out_dir,'ShearDesign.fig'));
fprintf('\nFigure saved to: %s/ShearDesign.png\n', out_dir);

fprintf('\n====================================================\n');
fprintf('Shear design complete.\n');
fprintf('====================================================\n');
