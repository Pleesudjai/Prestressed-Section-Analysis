function generateReport(beam, section, materials, prestress, loads, results_all, output_dir)
% GENERATEREPORT  Step-by-step calculation report saved as a .txt file.
%   Called from main_PrestressedBeamAnalysis after all stages complete.
%
%   FORMAT:  Professional exam-quality, symbolic equation then plugged-in values.
%   OUTPUT:  output_dir/Calculation_Report.txt

if nargin < 7 || isempty(output_dir);  output_dir = 'output';  end
if ~exist(output_dir, 'dir');  mkdir(output_dir);  end

rpt = fullfile(output_dir, 'Calculation_Report.txt');
fid = fopen(rpt, 'w');
if fid < 0;  error('Cannot open: %s', rpt);  end

%% ── helpers ──────────────────────────────────────────────────────────────
SEP  = @()     fprintf(fid, '%s\n', repmat('=',1,66));
sep  = @()     fprintf(fid, '%s\n', repmat('-',1,66));
BR   = @()     fprintf(fid, '\n');
PRT  = @(s)    fprintf(fid, '%s\n', s);
PRTf = @(fmt, varargin) fprintf(fid, [fmt '\n'], varargin{:});
PASS = @(v, lo, hi) ternary(v >= lo && v <= hi, 'PASS', ...
                             sprintf('FAIL  (limit [%.3f, +%.3f])', lo, hi));

%% ── derived quantities ───────────────────────────────────────────────────
A   = section.A;    Ic  = section.Ix;
yc  = section.yc;   yt  = section.yt;   yb  = section.yb;
St  = Ic / yt;      Sb  = Ic / yb;
r2  = Ic / A;       kb  = r2 / yt;      kt  = r2 / yb;
L   = beam.L;

% Self-weight (auto-computed inside analyzePrestressedBeam; pull from results)
w_sw  = results_all{1}.loads.self_weight;   % kip/in
w_SDL = 0;
if ~isempty(loads.distributed)
    w_SDL = abs(loads.distributed(1,3));
end
w_LL = 0;
if ~isempty(loads.distributed_live)
    w_LL = abs(loads.distributed_live(1,3));
end

% Prestress
n_tend    = length(prestress.tendons);
Aps_total = sum(cellfun(@(t) t.Aps, prestress.tendons));
Fi_total  = sum(cellfun(@(t) t.Aps * t.fpi, prestress.tendons));
eta       = 1 - prestress.losses;
Fe_total  = eta * Fi_total;
fpi_val   = prestress.tendons{1}.fpi;
Aps_one   = prestress.tendons{1}.Aps;

fc  = materials.fc;   fci = materials.fci;

% x-plot locations
x_plot  = beam.x_plot_fractions * L;   % [0, 120, 240, 384] in
n_loc   = length(x_plot);

loc_label = @(k) sprintf('x = %.0f in  (%.1f ft, x/L = %.3f)', ...
    x_plot(k), x_plot(k)/12, x_plot(k)/L);

%% ════════════════════════════════════════════════════════════════════════
%%  HEADER
%% ════════════════════════════════════════════════════════════════════════
SEP();
PRT('   CEE 530 PRESTRESSED CONCRETE  —  PROJECT 1');
PRT('   Double-T Beam: Stress Analysis and Shear Design');
SEP();
PRTf('   Student      : Chidchanok Pleesudjai (Fen)');
PRTf('   Date         : %s', datestr(now,'yyyy-mm-dd'));
PRTf('   Code edition : %s', materials.code_edition);
SEP();  BR();

%% ════════════════════════════════════════════════════════════════════════
%%  PART 1 — GIVEN DATA
%% ════════════════════════════════════════════════════════════════════════
SEP();  PRT('   PART 1 — GIVEN DATA');  SEP();  BR();

PRT('1.1  BEAM GEOMETRY');  sep();
PRTf('  Span:  L = (20 + 24 + 20) ft x 12 in/ft = %.0f in  (%.1f ft)', L, L/12);
BR();

PRT('1.2  CROSS-SECTION  (Double-T / TT, symmetric)');  sep();
PRT('  Section origin : y = 0 at BOTTOM of stems, positive upward');
PRT('  Top flange     : 120 in wide x 2 in thick   (y = 26 to 28 in)');
PRT('  Left stem      : x = -32 to -28.25 in        (y = 0 to 26 in)');
PRT('  Right stem     : x = +28.25 to +32 in        (y = 0 to 26 in)');
PRTf('  Total depth    : h = %.0f in', ...
    max(section.vertices(:,2))-min(section.vertices(:,2)));
BR();

PRT('1.3  MATERIAL PROPERTIES');  sep();
PRTf('  Concrete : f''c  = %.1f ksi   (28-day strength)', fc);
PRTf('             f''ci = %.1f ksi   (strength at transfer)', fci);
PRTf('             Ec   = %.0f ksi', materials.Ec);
PRTf('  PS steel : fpu  = %.0f ksi   (ultimate)', materials.fpu);
PRTf('             fpy  = %.0f ksi   (0.90 fpu)', materials.fpy);
PRTf('             Eps  = %.0f ksi', materials.Eps);
PRTf('  Mild     : fy   = %.0f ksi,  Es = %.0f ksi', materials.fy, materials.Es);
BR();

PRT('1.4  APPLIED LOADS');  sep();
PRT('  Self-weight:');
PRTf('    w_sw = Ac * rho_c');
PRTf('         = %.4f in^2 * (150/1728 lb/in^3) / 1000  [kip conversion]', A);
PRTf('         = %.6f kip/in  =  %.4f kip/ft', w_sw, w_sw*12);
BR();
PRT('  Superimposed Dead Load  (2-in concrete topping, 150 pcf):');
PRTf('    w_SDL = 2 in * 120 in * (150/1728 lb/in^3) / 1000');
PRTf('          = %.6f kip/in  =  %.4f kip/ft', w_SDL, w_SDL*12);
BR();
PRT('  Live Load:');
PRTf('    w_LL  = 0.42 kip/ft / 12 in/ft');
PRTf('          = %.6f kip/in  =  %.4f kip/ft', w_LL, w_LL*12);
BR();

PRT('1.5  PRESTRESSING TENDONS  (4 strands, fully bonded)');  sep();
PRTf('  Strand area    : Aps = %.3f in^2  per strand  (0.5-in diameter, Grade 270)', Aps_one);
PRTf('  Initial force  : Pi  = %.1f kips  per strand', Aps_one * fpi_val);
PRTf('    fpi = Pi/Aps = %.1f / %.3f = %.1f ksi', Aps_one*fpi_val, Aps_one, fpi_val);
BR();
PRT('  Tendon layout (y from bottom):');
PRT('    Tendon 1 — left stem,  straight : y = 6 in  (constant)');
PRT('    Tendon 2 — left stem,  harped   : y = 20.36 in at supports,');
PRT('                                       y = 6 in  from x = 240 to 528 in');
PRT('    Tendon 3 — right stem, straight : y = 6 in  (constant, mirror of T1)');
PRT('    Tendon 4 — right stem, harped   : same profile as Tendon 2');
BR();
PRTf('  Total initial prestress:   Fi = %d x %.1f = %.1f kips', ...
    n_tend, Aps_one*fpi_val, Fi_total);
PRTf('  Prestress losses:          %.0f%%  (given)', prestress.losses*100);
PRTf('  Loss factor:               eta = 1 - %.2f = %.2f', prestress.losses, eta);
PRTf('  Effective prestress:       Fe  = eta * Fi = %.2f * %.1f = %.1f kips', ...
    eta, Fi_total, Fe_total);
BR();

%% ════════════════════════════════════════════════════════════════════════
%%  PART 2 — SECTION PROPERTIES
%% ════════════════════════════════════════════════════════════════════════
SEP();  PRT('   PART 2 — SECTION PROPERTIES  (Shoelace / Green''s Theorem)');  SEP();  BR();

PRT('2.1  CROSS-SECTIONAL AREA');  sep();
PRT('  Formula (Shoelace):');
PRT('    A = (1/2) * |Sum_{i=1}^{n} (x_i * y_{i+1} - x_{i+1} * y_i)|');
BR();
PRTf('  A = %.4f in^2', A);
BR();

PRT('2.2  CENTROID  (from bottom, y = 0 at stem soffit)');  sep();
PRT('  Formula:');
PRT('    yc = [Sum (x_i*y_{i+1}-x_{i+1}*y_i)*(y_i+y_{i+1})] / (6*A)');
BR();
PRTf('  yc = %.4f in   from bottom', yc);
PRTf('  yt = y_top - yc  =  28.000 - %.4f  =  %.4f in   (centroid to TOP fiber)', yc, yt);
PRTf('  yb = yc - 0      =  %.4f - 0      =  %.4f in   (centroid to BOTTOM fiber)', yc, yb);
BR();

PRT('2.3  MOMENT OF INERTIA  (centroidal x-axis)');  sep();
PRT('  Formula:');
PRT('    Ic = (1/12)*|Sum (x_i*y_{i+1}-x_{i+1}*y_i)*(y_i^2+y_i*y_{i+1}+y_{i+1}^2)|');
BR();
PRTf('  Ic = %.2f in^4', Ic);
BR();

PRT('2.4  SECTION MODULI');  sep();
PRTf('  St = Ic / yt  =  %.2f / %.4f  =  %.4f in^3   (top)', Ic, yt, St);
PRTf('  Sb = Ic / yb  =  %.2f / %.4f  =  %.4f in^3   (bottom)', Ic, yb, Sb);
BR();

PRT('2.5  KERN POINTS  (limits for no-tension design)');  sep();
PRTf('  r^2 = Ic / A  =  %.2f / %.4f  =  %.4f in^2', Ic, A, r2);
PRTf('  kb  = r^2 / yt  =  %.4f / %.4f  =  %.4f in   (upper kern, bottom tendons)', r2, yt, kb);
PRTf('  kt  = r^2 / yb  =  %.4f / %.4f  =  %.4f in   (lower kern, top tendons)', r2, yb, kt);
BR();

%% ════════════════════════════════════════════════════════════════════════
%%  PART 3 — ALLOWABLE STRESSES
%% ════════════════════════════════════════════════════════════════════════
SEP();  PRTf('   PART 3 — ALLOWABLE STRESSES  [%s]', materials.code_edition);  SEP();  BR();

PRT('3.1  AT TRANSFER  (f''ci = 4.8 ksi)');  sep();
PRTf('  Compression (general):');
PRTf('    f_ci_allow = +0.60 * f''ci  =  0.60 * %.1f  =  +%.3f ksi', fci, materials.f_ci_allow);
PRTf('  Compression (end regions):');
if strcmp(materials.code_edition, 'ACI-318-19')
    PRTf('    f_ci_end   = +0.70 * f''ci  =  0.70 * %.1f  =  +%.3f ksi   [ACI-318-19]', ...
        fci, materials.f_ci_allow_end);
else
    PRTf('    f_ci_end   = +0.60 * f''ci  =  0.60 * %.1f  =  +%.3f ksi   [CEE530: same as general]', ...
        fci, materials.f_ci_allow_end);
end
PRTf('  Tension (general):');
PRTf('    f_ti_allow = -3 * sqrt(f''ci_psi)  =  -3 * sqrt(%.0f)  =  %.3f ksi', ...
    fci*1000, materials.f_ti_allow);
PRTf('  Tension (end regions):');
PRTf('    f_ti_end   = -6 * sqrt(f''ci_psi)  =  -6 * sqrt(%.0f)  =  %.3f ksi', ...
    fci*1000, materials.f_ti_allow_end);
BR();

PRT('3.2  AT SERVICE  (f''c = 6.0 ksi)');  sep();
PRTf('  Compression (sustained SW+SDL):');
PRTf('    f_cs_sust  = +0.45 * f''c  =  0.45 * %.1f  =  +%.3f ksi', fc, materials.f_cs_allow_sust);
PRTf('  Compression (total SW+SDL+LL):');
PRTf('    f_cs_total = +0.60 * f''c  =  0.60 * %.1f  =  +%.3f ksi', fc, materials.f_cs_allow_total);
PRTf('  Tension Class C:');
PRTf('    f_ts = -12 * sqrt(f''c_psi)  =  -12 * sqrt(%.0f)  =  %.3f ksi', ...
    fc*1000, materials.f_ts_allow);
PRTf('  Tension Class U boundary:');
if strcmp(materials.code_edition, 'ACI-318-19')
    PRTf('    f_tu = -7.5 * sqrt(f''c_psi)  =  -7.5 * sqrt(%.0f)  =  %.3f ksi   [ACI-318-19]', ...
        fc*1000, materials.f_tu_allow);
else
    PRTf('    f_tu = -6.0 * sqrt(f''c_psi)  =  -6.0 * sqrt(%.0f)  =  %.3f ksi   [CEE530]', ...
        fc*1000, materials.f_tu_allow);
end
BR();

%% ════════════════════════════════════════════════════════════════════════
%%  PART 4 — TENDON ECCENTRICITY
%% ════════════════════════════════════════════════════════════════════════
SEP();  PRT('   PART 4 — TENDON ECCENTRICITY AT ANALYSIS SECTIONS');  SEP();  BR();
PRT('  e = yc - y_tendon   (positive = tendon BELOW centroid)');
PRTf('  yc = %.4f in  (from bottom)', yc);
BR();

for k = 1:n_loc
    x_k = x_plot(k);
    [~, ik] = min(abs(beam.x - x_k));
    PRTf('  %-40s', loc_label(k));
    e_sum = 0;  P_sum = 0;
    for tt = 1:n_tend
        t     = prestress.tendons{tt};
        y_t   = section.yc - t.e(ik);
        e_t   = t.e(ik);
        PRTf('    Tendon %d: y_t = %.2f in  -->  e = yc - y_t = %.4f - %.2f = %+.4f in', ...
            tt, y_t, yc, y_t, e_t);
        e_sum = e_sum + (t.Aps * t.fpi) * e_t;
        P_sum = P_sum + t.Aps * t.fpi;
    end
    e_eff = e_sum / P_sum;
    PRTf('    Force-weighted average:  e_eff = %.4f in', e_eff);
    BR();
end

%% ════════════════════════════════════════════════════════════════════════
%%  PARTS 5-7 — STRESS ANALYSIS BY STAGE
%% ════════════════════════════════════════════════════════════════════════
part_num = 5;
for si = 1:length(results_all)
    ri       = results_all{si};
    sname    = ri.stage_name;
    f_C      = ri.stresses.fc_allow_compression;
    f_T      = ri.stresses.fc_allow_tension;
    mid_idx  = round(length(ri.x)/2);
    eta_stg  = ri.P(mid_idx) / Fi_total;

    SEP();
    switch sname
        case 'Transfer'
            PRTf('   PART %d — STAGE 1: TRANSFER  (eta = 1.00)', part_num);
            PRTf('   Concrete strength = f''ci = %.1f ksi  |  M = M_sw only', fci);
        case 'Service_Sustained'
            PRTf('   PART %d — STAGE 2: SERVICE - SUSTAINED  (eta = %.2f)', part_num, eta);
            PRTf('   Concrete strength = f''c = %.1f ksi  |  M = M_sw + M_SDL', fc);
        case 'Service_Total'
            PRTf('   PART %d — STAGE 3: SERVICE - TOTAL  (eta = %.2f)', part_num, eta);
            PRTf('   Concrete strength = f''c = %.1f ksi  |  M = M_sw + M_SDL + M_LL', fc);
    end
    SEP();  BR();

    % --- Equations ---
    PRTf('  SIGN CONVENTION:  Compression = positive (+)   Tension = negative (-)');
    BR();
    PRTf('  FIBER STRESS EQUATIONS:');
    BR();
    PRTf('    f_top = +F/Ac  -  F*e/St  +  M/St');
    PRTf('    f_bot = +F/Ac  +  F*e/Sb  -  M/Sb');
    BR();
    PRTf('  TERMS:');
    PRTf('    F   = effective prestress force at section  (kips)');
    PRTf('    Ac  = %.4f in^2', A);
    PRTf('    e   = yc - y_tendon(eff)  (in, positive when tendon below centroid)');
    PRTf('    St  = Ic/yt  =  %.2f / %.4f  =  %.4f in^3', Ic, yt, St);
    PRTf('    Sb  = Ic/yb  =  %.2f / %.4f  =  %.4f in^3', Ic, yb, Sb);
    PRTf('    M   = total moment at section  (kip-in, sagging = positive)');
    BR();

    if strcmp(sname, 'Transfer')
        PRTf('  F = Fi = %.1f kips  (full prestress, eta = 1.00, no losses at transfer)', Fi_total);
    else
        PRTf('  F = eta * Fi  =  %.2f * %.1f  =  %.1f kips', eta, Fi_total, Fe_total);
    end
    PRTf('  Allowable:  %.3f ksi (tension)  <=  f  <=  +%.3f ksi (compression)', f_T, f_C);
    BR();

    % --- Per-location ---
    for k = 1:n_loc
        x_k   = x_plot(k);
        [~,ik] = min(abs(ri.x - x_k));

        P_k  = ri.P(ik);
        e_k  = ri.e(ik);
        M_k  = ri.M(ik);
        ftp  = ri.stresses.f_top_prestress(ik);   % +P/A - P*e/St
        flp  = ri.stresses.f_top(ik);             % +M/St
        ftt  = ri.stresses.f_top_total(ik);
        fbp  = ri.stresses.f_bot_prestress(ik);   % +P/A + P*e/Sb
        fbl  = ri.stresses.f_bot(ik);             % -M/Sb
        fbt  = ri.stresses.f_bot_total(ik);

        % Individual numeric terms for display
        tPA  =  P_k / A;
        tTs  = -P_k * e_k / St;    % eccentricity term, top
        tTM  =  M_k / St;          % moment term, top
        tBs  =  P_k * e_k / Sb;    % eccentricity term, bottom
        tBM  = -M_k / Sb;          % moment term, bottom

        sep();
        PRTf('  SECTION:  %s', loc_label(k));
        sep();
        BR();

        % Moment
        switch sname
            case 'Transfer'
                PRTf('  Self-weight moment:');
                PRTf('    M_sw = w_sw * x * (L - x) / 2');
                PRTf('         = %.6f * %.1f * (%.0f - %.1f) / 2', w_sw, x_k, L, x_k);
                PRTf('         = %.2f kip-in  =  %.3f kip-ft', M_k, M_k/12);
            case 'Service_Sustained'
                PRTf('  Sustained moment  (M_sw + M_SDL):');
                PRTf('    M = (w_sw + w_SDL) * x * (L - x) / 2');
                PRTf('      = (%.6f + %.6f) * %.1f * (%.0f - %.1f) / 2', w_sw, w_SDL, x_k, L, x_k);
                PRTf('      = (%.6f) * %.1f * %.1f / 2', w_sw+w_SDL, x_k, L-x_k);
                PRTf('      = %.2f kip-in  =  %.3f kip-ft', M_k, M_k/12);
            case 'Service_Total'
                PRTf('  Total moment  (M_sw + M_SDL + M_LL):');
                PRTf('    M = (w_sw + w_SDL + w_LL) * x * (L - x) / 2');
                PRTf('      = (%.6f + %.6f + %.6f) * %.1f * (%.0f - %.1f) / 2', ...
                    w_sw, w_SDL, w_LL, x_k, L, x_k);
                PRTf('      = (%.6f) * %.1f * %.1f / 2', w_sw+w_SDL+w_LL, x_k, L-x_k);
                PRTf('      = %.2f kip-in  =  %.3f kip-ft', M_k, M_k/12);
        end
        BR();

        PRTf('  Prestress at this section:');
        PRTf('    F   = %.4f kips', P_k);
        PRTf('    e   = %.4f in   (force-weighted average eccentricity)', e_k);
        BR();

        % TOP FIBER
        PRT('  *** TOP FIBER ***');
        BR();
        PRT('    f_top  =  +F/Ac  -  F*e/St  +  M/St');
        BR();
        PRTf('           =  +%.4f/%.4f  -  %.4f*%.4f/%.4f  +  %.4f/%.4f', ...
            P_k, A, P_k, e_k, St, M_k, St);
        BR();
        PRTf('           =  +%.4f        -  (%.4f)             +  (%.4f)', ...
            tPA, abs(tTs), tTM);
        PRTf('              [P/A unif]      [eccen, top]           [moment]');
        BR();
        PRTf('    f_top  = %+.4f ksi', ftt);
        BR();
        PRTf('    Allowable:  %.3f ksi (T)  <=  f_top  <=  +%.3f ksi (C)', f_T, f_C);
        PRTf('    RESULT:  f_top = %+.4f ksi   -->  %s', ftt, PASS(ftt, f_T, f_C));
        BR();

        % BOTTOM FIBER
        PRT('  *** BOTTOM FIBER ***');
        BR();
        PRT('    f_bot  =  +F/Ac  +  F*e/Sb  -  M/Sb');
        BR();
        PRTf('           =  +%.4f/%.4f  +  %.4f*%.4f/%.4f  -  %.4f/%.4f', ...
            P_k, A, P_k, e_k, Sb, M_k, Sb);
        BR();
        PRTf('           =  +%.4f        +  (%.4f)             -  (%.4f)', ...
            tPA, tBs, abs(tBM));
        PRTf('              [P/A unif]      [eccen, bot]           [moment]');
        BR();
        PRTf('    f_bot  = %+.4f ksi', fbt);
        BR();
        PRTf('    Allowable:  %.3f ksi (T)  <=  f_bot  <=  +%.3f ksi (C)', f_T, f_C);
        PRTf('    RESULT:  f_bot = %+.4f ksi   -->  %s', fbt, PASS(fbt, f_T, f_C));
        BR();
    end

    part_num = part_num + 1;
end

%% ════════════════════════════════════════════════════════════════════════
%%  SHEAR DESIGN PART
%% ════════════════════════════════════════════════════════════════════════
SEP();  PRTf('   PART %d — SHEAR DESIGN  [ACI 318-19, Naaman Sec. 5.7]', part_num);  SEP();  BR();

lambda = 1.0;   phi_v = 0.75;
Av_bar = 2 * 0.11;   fy_s = 60.0;
sqrtfc = sqrt(fc * 1000) / 1000;

h_sec   = max(section.vertices(:,2)) - min(section.vertices(:,2));
y_minv  = min(section.vertices(:,2));
xv      = sort(section.vertices(abs(section.vertices(:,2)-y_minv)<1e-6, 1));
bw = 0;
for kv = 1:2:length(xv)-1;  bw = bw + (xv(kv+1)-xv(kv));  end

xg = beam.x;  ng = length(xg);
y_ps = zeros(1,ng);
for i = 1:n_tend
    t = prestress.tendons{i};
    y_ps = y_ps + t.Aps * t.y;
end
y_ps = y_ps / Aps_total;
dp   = max(h_sec - y_ps, 0.8*h_sec);
e_ps = yc - y_ps;
Pe   = Fe_total * ones(1,ng);

Vp = zeros(1,ng);
for i = 1:n_tend
    t = prestress.tendons{i};
    if isfield(t,'y_profile') && ~isempty(t.y_profile)
        dy_dx = zeros(1,ng);
        xp = t.x_profile;  yp = t.y_profile;
        for seg = 1:length(xp)-1
            msk = (xg>=xp(seg)) & (xg<=xp(seg+1));
            dy_dx(msk) = (yp(seg+1)-yp(seg))/(xp(seg+1)-xp(seg));
        end
        Pe_i = t.Aps * t.fpi * eta;
        Vp   = Vp + Pe_i * sin(atan(abs(dy_dx)));
    end
end

w_DL = w_sw + w_SDL;
w_u  = 1.2*w_DL + 1.6*w_LL;
w_ext= 1.2*w_SDL + 1.6*w_LL;

Vd   = w_DL * (L/2 - xg);
Md   = w_DL * xg .* (L-xg) / 2;
Vu   = abs(w_u * (L/2 - xg));
Vi   = abs(w_ext * (L/2 - xg));
Mmax = w_ext * xg .* (L-xg) / 2;
Mmax(Mmax<1e-6) = 1e-6;

fce    = Pe/A + Pe .* e_ps * yb ./ Ic;
fd_bot = abs(Md) * yb ./ Ic;
fr_sh  = 6 * lambda * sqrtfc;
Mcr    = max((Ic/yb)*(fr_sh + fce - fd_bot), 0);

Vci_t1 = 0.6*lambda*sqrtfc*bw*dp;
Vci_t2 = abs(Vd);
Vci_t3 = Vi.*Mcr./Mmax;  Vci_t3(~isfinite(Vci_t3))=1e9;
Vci    = max(Vci_t1+Vci_t2+Vci_t3, 1.7*lambda*sqrtfc*bw*dp);
Vcw    = (3.5*lambda*sqrtfc + 0.3*Pe/A).*bw.*dp + Vp;
Vc     = min(Vci, Vcw);

Vs_req = max(Vu/phi_v - Vc, 0);
Avs_d  = Vs_req ./ (fy_s * dp);
Avs_c1 = 0.75*sqrtfc*bw/fy_s;
Avs_c2 = 50/1000*bw/fy_s;
Avs_c3 = Aps_total*materials.fpu./(80*fy_s*dp).*sqrt(dp/bw);
Avs_min= max([Avs_c1, Avs_c2, max(Avs_c3)]);
Avs    = max(Avs_d, Avs_min);

s_lim  = min(3*h_sec/4, 24);
s_max  = s_lim * ones(1,ng);
s_req  = Av_bar ./ Avs;
s_des  = max(floor(min(s_req,s_max)),1);

[~,ic] = min(abs(xg - max(dp(1), 0.8*h_sec)));
mgd    = round(ng/2);

%% ── Shear equations ──
PRT('  EQUATIONS:');  sep();
PRT('    Vci = 0.6*lambda*sqrt(f''c)*bw*dp + Vd + (Vi/Mmax)*Mcr');
PRT('          >= Vci_min = 1.7*lambda*sqrt(f''c)*bw*dp');
PRT('    Mcr = (Ic/yb)*(fr + fce - fd)');
PRTf('    fr  = 6*lambda*sqrt(f''c_psi)/1000');
PRT('    fce = Pe/Ac + Pe*e*yb/Ic   (at bottom fiber)');
PRT('    fd  = Md*yb/Ic             (DL stress at bottom, tensile magnitude)');
BR();
PRT('    Vcw = (3.5*lambda*sqrt(f''c) + 0.3*fpc)*bw*dp + Vp');
PRT('    fpc = Pe/Ac                (at centroid)');
PRT('    Vp  = sum Pe_i * sin(theta_i)  (vertical PS component, upward)');
BR();
PRT('    Vc  = min(Vci, Vcw)');
PRT('    Vs_req = max(Vu/phi - Vc, 0)');
PRT('    Av/s = Vs/(fy*dp)   for vertical stirrups (theta_s = 90 deg)');
BR();

%% ── Section properties for shear ──
PRT('  SECTION PROPERTIES FOR SHEAR:');  sep();
PRTf('    h   = %.2f in', h_sec);
PRTf('    bw  = %.2f in   (sum of stem widths at bottom)', bw);
PRTf('    Ac  = %.4f in^2', A);
PRTf('    Ic  = %.2f in^4', Ic);
PRTf('    yb  = %.4f in', yb);
PRTf('    yc  = %.4f in   (from bottom)', yc);
BR();

%% ── Factored loads ──
PRT('  FACTORED LOADS:');  sep();
PRTf('    wu = 1.2*wDL + 1.6*wLL');
PRTf('       = 1.2*(w_sw + w_SDL) + 1.6*w_LL');
PRTf('       = 1.2*(%.6f + %.6f) + 1.6*(%.6f)', w_sw, w_SDL, w_LL);
PRTf('       = 1.2*(%.6f) + 1.6*(%.6f)', w_DL, w_LL);
PRTf('       = %.6f kip/in  =  %.4f kip/ft', w_u, w_u*12);
BR();
PRTf('    Max Vu = wu*L/2 = %.6f * %.0f/2 = %.3f kips', w_u, L, max(Vu));
PRTf('    Max Mu = wu*L^2/8 = %.6f * %.0f^2/8 = %.2f kip-in = %.2f kip-ft', ...
    w_u, L, max(Vu.*(L/2-xg)+Vc*0), w_u*L^2/8/12);
BR();

%% ── Critical section ──
PRTf('  CRITICAL SECTION:  x = %.2f in = %.3f ft  (dp from face of support)', ...
    xg(ic), xg(ic)/12);
PRTf('    dp = max(h-y_ps, 0.8h) = max(%.3f, %.3f) = %.3f in', ...
    h_sec-y_ps(ic), 0.8*h_sec, dp(ic));
BR();

%% ── Step 1: fr ──
PRT('  STEP 1 — MODULUS OF RUPTURE FOR SHEAR  (fr)');  sep();
PRTf('    fr = 6 * lambda * sqrt(f''c_psi) / 1000');
PRTf('       = 6 * %.1f * sqrt(%.0f) / 1000', lambda, fc*1000);
PRTf('       = 6 * %.1f * %.4f', lambda, sqrtfc);
PRTf('       = %.4f ksi', fr_sh);
BR();

%% ── Step 2: fce ──
PRT('  STEP 2 — PRESTRESS STRESS AT BOTTOM FIBER  (fce)');  sep();
PRT('    fce = Pe/Ac  +  Pe*e*yb/Ic');
BR();
PRTf('        = %.3f/%.4f  +  %.3f*%.4f*%.4f/%.2f', ...
    Fe_total, A, Fe_total, e_ps(ic), yb, Ic);
PRTf('        = %.4f  +  %.4f', Fe_total/A, Fe_total*e_ps(ic)*yb/Ic);
PRTf('        = %.4f ksi   (compression, +)', fce(ic));
BR();

%% ── Step 3: fd ──
PRT('  STEP 3 — DEAD LOAD FLEXURE STRESS AT BOTTOM FIBER  (fd)');  sep();
PRT('    fd = Md * yb / Ic   (tensile magnitude, used in Mcr formula)');
BR();
PRTf('    Md = w_DL * x * (L-x) / 2');
PRTf('       = %.6f * %.2f * (%.0f - %.2f) / 2', w_DL, xg(ic), L, xg(ic));
PRTf('       = %.3f kip-in', Md(ic));
BR();
PRTf('    fd = %.3f * %.4f / %.2f', Md(ic), yb, Ic);
PRTf('       = %.4f ksi   (tensile magnitude)', fd_bot(ic));
BR();

%% ── Step 4: Mcr ──
PRT('  STEP 4 — CRACKING MOMENT  (Mcr)');  sep();
PRT('    Mcr = (Ic/yb) * (fr + fce - fd)');
BR();
PRTf('        = (%.2f / %.4f) * (%.4f + %.4f - %.4f)', Ic, yb, fr_sh, fce(ic), fd_bot(ic));
PRTf('        = %.4f * %.4f', Ic/yb, fr_sh+fce(ic)-fd_bot(ic));
PRTf('        = %.2f kip-in  =  %.2f kip-ft', Mcr(ic), Mcr(ic)/12);
if (Ic/yb)*(fr_sh+fce(ic)-fd_bot(ic)) < 0
    PRT('    NOTE: Raw Mcr < 0  -->  set Mcr = 0  (fd exceeds fr + fce)');
    PRT('          Vci_min governs at this location.');
end
BR();

%% ── Step 5: Vci ──
PRT('  STEP 5 — FLEXURAL-SHEAR STRENGTH  (Vci)');  sep();
PRT('    Vci = 0.6*lambda*sqrt(f''c)*bw*dp  +  Vd  +  (Vi/Mmax)*Mcr');
BR();
PRTf('    w_ext = 1.2*w_SDL + 1.6*w_LL');
PRTf('          = 1.2*%.6f + 1.6*%.6f  =  %.6f kip/in', w_SDL, w_LL, w_ext);
PRTf('    Vi    = w_ext*(L/2 - x)  =  %.6f*(%.1f-%.2f)  =  %.3f kips', ...
    w_ext, L/2, xg(ic), Vi(ic));
PRTf('    Mmax  = w_ext*x*(L-x)/2  =  %.6f*%.2f*(%.0f-%.2f)/2  =  %.2f kip-in', ...
    w_ext, xg(ic), L, xg(ic), Mmax(ic));
PRTf('    Vd    = w_DL*(L/2-x)     =  %.6f*(%.1f-%.2f)  =  %.3f kips', ...
    w_DL, L/2, xg(ic), Vci_t2(ic));
BR();
PRTf('    Term 1: 0.6 * %.1f * %.4f * %.2f * %.3f  =  %.3f kips', ...
    lambda, sqrtfc, bw, dp(ic), Vci_t1(ic));
PRTf('    Term 2: Vd                               =  %.3f kips', Vci_t2(ic));
PRTf('    Term 3: (%.3f/%.2f) * %.2f             =  %.3f kips', ...
    Vi(ic), Mmax(ic), Mcr(ic), Vci_t3(ic));
BR();
PRTf('    Vci (sum)   = %.3f + %.3f + %.3f  =  %.3f kips', ...
    Vci_t1(ic), Vci_t2(ic), min(Vci_t3(ic),999), Vci_t1(ic)+Vci_t2(ic)+min(Vci_t3(ic),999));
PRTf('    Vci_min = 1.7*%.1f*%.4f*%.2f*%.3f   =  %.3f kips', ...
    lambda, sqrtfc, bw, dp(ic), 1.7*lambda*sqrtfc*bw*dp(ic));
PRTf('    Vci (governing) = max(above)          =  %.3f kips', Vci(ic));
BR();

%% ── Step 6: Vcw ──
PRT('  STEP 6 — WEB-SHEAR STRENGTH  (Vcw)');  sep();
PRT('    Vcw = (3.5*lambda*sqrt(f''c) + 0.3*fpc) * bw * dp  +  Vp');
BR();
PRTf('    fpc = Pe/Ac  =  %.3f/%.4f  =  %.4f ksi   (at centroid)', Fe_total, A, Fe_total/A);
PRTf('    Vp  = %.3f kips   (upward component, harped tendons at x = %.2f in)', Vp(ic), xg(ic));
BR();
PRTf('    (3.5*%.1f*%.4f + 0.3*%.4f) * %.2f * %.3f', ...
    lambda, sqrtfc, Fe_total/A, bw, dp(ic));
PRTf('    = (%.4f + %.4f) * %.2f * %.3f', 3.5*lambda*sqrtfc, 0.3*Fe_total/A, bw, dp(ic));
PRTf('    = %.4f * %.2f * %.3f', 3.5*lambda*sqrtfc+0.3*Fe_total/A, bw, dp(ic));
PRTf('    = %.3f kips', (3.5*lambda*sqrtfc+0.3*Fe_total/A)*bw*dp(ic));
BR();
PRTf('    Vcw = %.3f + %.3f  =  %.3f kips', ...
    (3.5*lambda*sqrtfc+0.3*Fe_total/A)*bw*dp(ic), Vp(ic), Vcw(ic));
BR();

%% ── Step 7: Vc ──
PRT('  STEP 7 — CONCRETE SHEAR STRENGTH  Vc = min(Vci, Vcw)');  sep();
PRTf('    Vc  = min(%.3f, %.3f)  =  %.3f kips', Vci(ic), Vcw(ic), Vc(ic));
PRTf('    phi*Vc  =  %.2f * %.3f  =  %.3f kips', phi_v, Vc(ic), phi_v*Vc(ic));
BR();
PRTf('    Vu (at critical section)  =  %.3f kips', Vu(ic));
if Vu(ic) <= phi_v*Vc(ic)
    PRTf('    CHECK: Vu = %.3f <= phi*Vc = %.3f   -->  PASS (min stirrups only)', ...
        Vu(ic), phi_v*Vc(ic));
else
    PRTf('    CHECK: Vu = %.3f > phi*Vc = %.3f    -->  Stirrups REQUIRED', ...
        Vu(ic), phi_v*Vc(ic));
end
BR();

%% ── Step 8: Vs ──
PRT('  STEP 8 — REQUIRED STIRRUP CONTRIBUTION  (Vs)');  sep();
PRT('    Vs_req = max(Vu/phi - Vc, 0)');
PRTf('           = max(%.3f/%.2f - %.3f, 0)', Vu(ic), phi_v, Vc(ic));
PRTf('           = max(%.3f, 0)', Vu(ic)/phi_v - Vc(ic));
PRTf('           = %.3f kips', Vs_req(ic));
BR();

%% ── Step 9: Stirrups ──
PRTf('  STEP 9 — STIRRUP DESIGN  (#3 U-stirrups: Av = 2 x 0.11 = %.2f in^2)', Av_bar);  sep();
BR();
PRT('  Minimum Av/s  (ACI 318-19 Table 9.6.3.3, prestressed member):');
BR();
PRTf('    Criterion 1:  Av/s_min1 = 0.75*sqrt(f''c_psi)*bw/fy  (ksi system)');
PRTf('                           = 0.75 * %.4f * %.2f / %.0f', sqrtfc, bw, fy_s);
PRTf('                           = %.5f in^2/in', Avs_c1);
BR();
PRTf('    Criterion 2:  Av/s_min2 = (50 psi)*bw/fy  (psi->ksi: 50/1000)');
PRTf('                           = (50/1000) * %.2f / %.0f', bw, fy_s);
PRTf('                           = %.5f in^2/in', Avs_c2);
BR();
PRTf('    Criterion 3:  Av/s_min3 = Aps*fpu / (80*fy*dp) * sqrt(dp/bw)');
PRTf('                           = %.4f*%.0f / (80*%.0f*%.3f) * sqrt(%.3f/%.2f)', ...
    Aps_total, materials.fpu, fy_s, dp(ic), dp(ic), bw);
PRTf('                           = %.5f in^2/in', max(Avs_c3));
BR();
PRTf('    Av/s_min = max(%.5f, %.5f, %.5f)  =  %.5f in^2/in', ...
    Avs_c1, Avs_c2, max(Avs_c3), Avs_min);
BR();
PRT('  Required Av/s from Vs:');
PRTf('    Av/s_demand = Vs / (fy * dp)');
PRTf('               = %.3f / (%.0f * %.3f)', Vs_req(ic), fy_s, dp(ic));
PRTf('               = %.5f in^2/in', Avs_d(ic));
BR();
PRTf('    Av/s_governing = max(demand, min)  =  max(%.5f, %.5f)  =  %.5f in^2/in', ...
    Avs_d(ic), Avs_min, Avs(ic));
BR();
PRT('  Required spacing:');
PRTf('    s_req = Av / (Av/s)  =  %.2f / %.5f  =  %.2f in', Av_bar, Avs(ic), Av_bar/Avs(ic));
BR();
PRT('  Spacing limits:');
PRTf('    s_max (basic) = min(3h/4, 24 in)  =  min(%.2f, 24)  =  %.1f in', ...
    3*h_sec/4, s_lim);
BR();
PRTf('    s_design = min(s_req, s_max)  =  min(%.2f, %.1f)  =  %d in', ...
    Av_bar/Avs(ic), s_max(ic), s_des(ic));
PRTf('    USE:  #3 U-stirrups @ %d in  throughout region near support', s_des(ic));
BR();

%% ════════════════════════════════════════════════════════════════════════
%%  SUMMARY TABLE
%% ════════════════════════════════════════════════════════════════════════
part_num = part_num + 1;
SEP();  PRTf('   PART %d — STRESS SUMMARY  (all stages, all sections)', part_num);  SEP();  BR();

for si = 1:length(results_all)
    ri    = results_all{si};
    sname = ri.stage_name;
    f_C   = ri.stresses.fc_allow_compression;
    f_T   = ri.stresses.fc_allow_tension;

    PRTf('  Stage: %-20s   Limits: [%.3f, +%.3f] ksi', sname, f_T, f_C);
    PRTf('  %-24s  %+10s  %+10s  %-12s  %-12s', ...
        'Section', 'f_top (ksi)', 'f_bot (ksi)', 'Top status', 'Bot status');
    PRTf('  %s', repmat('-',1,72));
    for k = 1:n_loc
        x_k = x_plot(k);
        [~,ik] = min(abs(ri.x - x_k));
        ftt = ri.stresses.f_top_total(ik);
        fbt = ri.stresses.f_bot_total(ik);
        PRTf('  x = %5.0f in (%5.1f ft)   %+10.4f  %+10.4f  %-12s  %-12s', ...
            x_k, x_k/12, ftt, fbt, PASS(ftt,f_T,f_C), PASS(fbt,f_T,f_C));
    end
    BR();
end

%% ── End ──
SEP();  PRT('   END OF CALCULATION REPORT');  SEP();
fclose(fid);
fprintf('\n  Report saved to: %s\n', rpt);

end

%% ── local utility ──────────────────────────────────────────────────────
function s = ternary(cond, yes, no)
    if cond;  s = yes;  else;  s = no;  end
end
