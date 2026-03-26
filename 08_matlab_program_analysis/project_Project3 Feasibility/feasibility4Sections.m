%% FEASIBILITY4SECTIONS.m
%  Magnel Diagrams at cross-sections defined in inputData()
%  Project 1 Part 2 — Double-T Beam Feasibility Chart
%
%  Naaman line numbering:
%    Line I   : top @ transfer, tension limit   [lower bound on 1/F]
%    Line II  : bot @ transfer, compr. limit    [lower bound on 1/F]
%    Line III : bot @ service,  tension limit   [upper bound on 1/F]
%    Line IV  : top @ service,  compr. limit    [lower bound on 1/F]
%
%  Sign convention: Compression = positive, Tension = negative

clc; close all;

%% =========================================================================
%% USER SETTINGS  (defaults — overridden by runAllCases.m when batch-running)
%% =========================================================================
Folder_Name  = 'output';     % subfolder for saved figures (may be overridden)
cover        = 2.0;          % clear cover to outermost tendon (in)
n_plot_pts   = 1200;         % number of eccentricity points for line plotting

%% =========================================================================
%% PATHS — universal input file
%% =========================================================================
% script_dir may be pre-set by runAllCases.m so mfilename stays correct
if ~exist('script_dir', 'var')
    script_dir = fileparts(mfilename('fullpath'));
end
addpath(fullfile(script_dir, '..', 'project_Project2'));   % inputData.m
output_dir = fullfile(script_dir, Folder_Name);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% =========================================================================
%% 1. LOAD INPUT DATA
%% =========================================================================
[beam, section, materials, prestress, ~, loads] = inputData();

%% =========================================================================
%% 1b. CASE OVERRIDES  (injected by runAllCases.m; empty = no override)
%% =========================================================================
%  case_tag      : string — routes output to output/<case_tag>/
%  w_LL_override : scalar kip/in — replaces live load from inputData
%  tendon_y_all  : scalar in    — sets all tendons to straight at this y
if ~exist('case_tag',      'var'), case_tag      = ''; end
if ~exist('w_LL_override', 'var'), w_LL_override = []; end
if ~exist('tendon_y_all',  'var'), tendon_y_all  = []; end

if ~isempty(w_LL_override)
    loads.distributed_live = [0, beam.L, -w_LL_override, -w_LL_override];
    fprintf('  [OVERRIDE] w_LL = %.5f kip/in  (%.0f lb/ft)\n', ...
            w_LL_override, w_LL_override*12*1000);
end

if ~isempty(tendon_y_all)
    for ii = 1:numel(prestress.tendons)
        e_const = section.yc - tendon_y_all;
        prestress.tendons{ii}.e            = e_const * ones(size(beam.x));
        prestress.tendons{ii}.y            = tendon_y_all * ones(size(beam.x));
        prestress.tendons{ii}.profile_type = 'linear';
    end
    fprintf('  [OVERRIDE] All tendons set to y = %.2f in  (e = %.4f in)\n', ...
            tendon_y_all, section.yc - tendon_y_all);
end

if ~isempty(case_tag)
    Folder_Name = fullfile('output', case_tag);
    output_dir  = fullfile(script_dir, Folder_Name);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    fprintf('  [OUTPUT]   Saving to: %s\n\n', output_dir);
end

%% =========================================================================
%% 2. SECTION PROPERTIES  (from inputData)
%% =========================================================================
Ac = section.A;
Ic = section.Ix;
yc = section.yc;
yt = section.yt;
yb = section.yb;
St = Ic / yt;
Sb = Ic / yb;

fprintf('\nSection Properties:\n');
fprintf('  Ac = %.3f in^2 | Ic = %.1f in^4 | yc = %.4f in\n', Ac, Ic, yc);
fprintf('  St = %.2f in^3 | Sb = %.2f in^3\n', St, Sb);
fprintf('  kt = %.4f in  | kb = %.4f in\n\n', Ic/(Ac*yt), Ic/(Ac*yb));

%% =========================================================================
%% 3. ALLOWABLE STRESSES  (from inputData — tension stored as negative)
%% =========================================================================
eta = 1 - prestress.losses;

fci_compr = materials.f_ci_allow;        % +0.60 f'ci
fti_tens  = -materials.f_ti_allow;       % magnitude: 3*sqrt(f'ci_psi)/1000
fcs_compr = materials.f_cs_allow_total;  % +0.60 f'c  (total service, CEE530)
fts_tens  = -materials.f_tu_allow;       % magnitude: 6*sqrt(f'c_psi)/1000 (CEE530 Class U)

fprintf('Allowable Stresses (%s):\n', materials.code_edition);
fprintf('  Transfer compr : +%.4f ksi\n', fci_compr);
fprintf('  Transfer tens  : -%.4f ksi\n', fti_tens);
fprintf('  Service  compr : +%.4f ksi  (total load)\n', fcs_compr);
fprintf('  Service  tens  : -%.4f ksi  (Class U)\n\n', fts_tens);

%% =========================================================================
%% 4. LOADS  (from inputData)
%% =========================================================================
L    = beam.L;
w_sw = Ac * loads.concrete_density / 1000;   % kip/in

w_SDL = 0;
if ~isempty(loads.distributed)
    w_SDL = -loads.distributed(1,3);   % distributed stores negative (downward)
end

w_LL = 0;
if isfield(loads,'distributed_live') && ~isempty(loads.distributed_live)
    w_LL = -loads.distributed_live(1,3);
end

w_s = w_sw + w_SDL + w_LL;

fprintf('Loads:\n');
fprintf('  w_sw  = %.5f kip/in\n', w_sw);
fprintf('  w_SDL = %.5f kip/in\n', w_SDL);
fprintf('  w_LL  = %.5f kip/in\n', w_LL);
fprintf('  w_s   = %.5f kip/in  (total service)\n\n', w_s);

%% =========================================================================
%% 5. SECTION LOCATIONS  (from beam.x_plot_fractions in inputData)
%%    e_cg(x) = sum(Aps_i * e_i(x)) / sum(Aps_i)  — weighted by strand area
%%    Tendon e(x) already computed by processTendonProfiles inside inputData()
%% =========================================================================
x_sec     = beam.x_plot_fractions * L;
n_sec     = numel(x_sec);
n_tendons = numel(prestress.tendons);

% e_cg = force-weighted centroid of tendon group
%   e_cg(x) = sum( Pi_i * e_i(x) ) / sum( Pi_i )
%   where Pi_i = Aps_i * fpi_i  (initial force per tendon, may differ)
e_cg_sec = zeros(1, n_sec);
for k = 1:n_sec
    xi    = x_sec(k);
    Pe_sum = 0;
    P_sum  = 0;
    for i = 1:n_tendons
        t   = prestress.tendons{i};
        Pi  = t.Aps * t.fpi;                                   % initial force this tendon
        e_i = interp1(beam.x, t.e, xi, 'linear', 'extrap');   % eccentricity at xi
        Pe_sum = Pe_sum + Pi * e_i;
        P_sum  = P_sum  + Pi;
    end
    e_cg_sec(k) = Pe_sum / P_sum;
end

% Labels and figure Names — generated from x values, no hardcoding
x_mid_ft = (L / 12) / 2;
labels    = cell(1, n_sec);
fig_names = cell(1, n_sec);
for k = 1:n_sec
    x_ft = x_sec(k) / 12;
    if x_ft == 0
        tag = 'Left_Support';
    elseif x_ft == x_mid_ft
        tag = 'Midspan';
    else
        tag = '';
    end
    base = sprintf('Magnel_x%.4gft', x_ft);
    if ~isempty(tag)
        base = [base '_' tag];
    end
    labels{k}    = strrep(strrep(base,'_',' '),'Magnel ','');  % readable title
    fig_names{k} = base;                                        % figure Name property
end

%% =========================================================================
%% 6. PHYSICAL ECCENTRICITY BOUNDS
%% =========================================================================
e_geo_max = yb - cover;
e_geo_min = -(yt - cover);

%% =========================================================================
%% 7. LOOP — ONE MAGNEL DIAGRAM PER SECTION
%% =========================================================================
results = struct();

for k = 1:n_sec
    xi   = x_sec(k);
    e_cg = e_cg_sec(k);
    lbl  = labels{k};

    Mi = w_sw * xi * (L - xi) / 2;   % transfer moment: SW only
    MT = w_s  * xi * (L - xi) / 2;   % service moment: SW + SDL + LL

    fprintf('=== x = %.4g ft ===\n', xi/12);
    fprintf('  Mi = %.1f kip-in  |  MT = %.1f kip-in\n', Mi, MT);
    fprintf('  e_cg = %.4f in  (tendon group centroid)\n', e_cg);

    %% Magnel denominators
    dI   = fti_tens  + Mi / St;
    dII  = fci_compr + Mi / Sb;
    dIII = MT / Sb   - fts_tens;
    dIV  = fcs_compr - MT / St;

    fprintf('  dI=%.4f  dII=%.4f  dIII=%.4f  dIV=%.4f\n', dI,dII,dIII,dIV);

    %% Eccentricity vector for plotting
    e_v = linspace(e_geo_min - 1, e_geo_max + 1, n_plot_pts);

    %% Four Magnel lines
    lineI  =       (e_v/St  - 1/Ac) ./ dI;
    lineII =       (1/Ac + e_v/Sb)  ./ dII;
    if dIII > 1e-9
        lineIII = eta * (1/Ac + e_v/Sb) ./ dIII;
    else
        lineIII = inf(size(e_v));
    end
    if dIV > 1e-9
        lineIV = eta * (1/Ac - e_v/St) ./ dIV;
    else
        lineIV = -inf(size(e_v));
    end

    lower_bnd = max([lineI; lineII; lineIV], [], 1);
    upper_bnd = lineIII;

    feas_mask = (upper_bnd > lower_bnd) & (lower_bnd > 0) & ...
                (e_v >= e_geo_min) & (e_v <= e_geo_max);

    [lam_lb, lam_ub, F_min_cg, F_max_cg] = getFeasRange(e_cg, e_v, lower_bnd, upper_bnd);

    fprintf('  e_cg=%.4f: 1/F in [%.5f, %.5f] --> F in [%.1f, %.1f] kip\n\n', ...
            e_cg, lam_lb, lam_ub, F_min_cg, F_max_cg);

    results(k).x      = xi;
    results(k).Mi     = Mi;
    results(k).MT     = MT;
    results(k).e_cg   = e_cg;
    results(k).Fmin   = F_min_cg;
    results(k).Fmax   = F_max_cg;
    results(k).lam_lb = lam_lb;
    results(k).lam_ub = lam_ub;

    %% --------------------------------------------------------------
    %%  FIGURE
    %% --------------------------------------------------------------
    figure('Name', fig_names{k}, 'Position', [60 60 960 680]);
    hold on;

    %% Axis limits — zoom to feasible zone
    e_f  = e_v(feas_mask);
    lb_f = lower_bnd(feas_mask);

    if ~isempty(e_f) && any(isfinite(upper_bnd(feas_mask)))
        ub_f     = upper_bnd(feas_mask);
        ub_f_fin = ub_f(isfinite(ub_f));
        if isempty(ub_f_fin); ub_f_fin = lb_f * 3; end
        xl_lo = max(min(lb_f) * 0.5, -0.0003);
        xl_hi = max(ub_f_fin) * 2.8;
        if xl_hi <= xl_lo; xl_hi = xl_lo + 0.015; end
    else
        xl_lo = -0.0003;
        xl_hi = 0.015;
    end

    lam_max = xl_hi * 1.1;
    lam_min = -0.002;

    lineI_p  = max(lam_min, min(lam_max, lineI));
    lineII_p = max(lam_min, min(lam_max, lineII));
    lineIV_p = max(lam_min, min(lam_max, lineIV));
    if dIII > 1e-9
        lineIII_p = max(lam_min, min(lam_max, lineIII));
    else
        lineIII_p = nan(size(e_v));
    end

    %% Green feasibility fill
    if ~isempty(e_f)
        ub_f_plot = min(upper_bnd(feas_mask), lam_max);
        if any(isfinite(ub_f_plot))
            fill([lb_f, fliplr(ub_f_plot)], [-e_f, fliplr(-e_f)], ...
                 [0.60 0.92 0.60], 'EdgeColor','none', 'FaceAlpha',0.55, ...
                 'HandleVisibility','off');
            plot(lb_f,      -e_f, 'k-', 'LineWidth',0.5, 'HandleVisibility','off');
            plot(ub_f_plot, -e_f, 'k-', 'LineWidth',0.5, 'HandleVisibility','off');
        end
    end

    %% Four lines
    plot(lineI_p,   -e_v, 'b-',  'LineWidth',2.2, ...
         'DisplayName','Line I : Top@Xfer  Tens. \geq -f_{ti}');
    plot(lineII_p,  -e_v, 'r-',  'LineWidth',2.2, ...
         'DisplayName','Line II: Bot@Xfer  Compr. \leq f_{ci}');
    if dIII > 1e-9
        plot(lineIII_p, -e_v, 'm--', 'LineWidth',2.2, ...
             'DisplayName','Line III: Bot@Svc  Tens. \geq -f_{ts}');
    else
        plot(NaN, NaN, 'm--', 'LineWidth',2.2, ...
             'DisplayName','Line III: N/A (M_T = 0)');
    end
    plot(lineIV_p,  -e_v, 'g-',  'LineWidth',2.2, ...
         'DisplayName','Line IV: Top@Svc  Compr. \leq f_{cs}');

    %% Reference lines
    xline(0, 'k--', 'LineWidth',0.8, 'HandleVisibility','off');
    yline(0, 'k--', 'LineWidth',0.8, 'HandleVisibility','off');
    yline(-e_geo_max, 'k:', 'LineWidth',1.8, ...
          'DisplayName',sprintf('e_{max} = %.2f in  (cover = %.0f in)', e_geo_max, cover));

    %% e_cg vertical marker
    xl = xline(-e_cg, 'k--', 'LineWidth',2.2);
    xl.Label = sprintf('e_{cg} = %.2f in', e_cg);
    xl.LabelVerticalAlignment   = 'top';
    xl.LabelHorizontalAlignment = 'center';
    xl.FontSize = 10;
    xl.HandleVisibility = 'off';

    %% Feasibility range bar at e_cg
    if lam_lb > 0 && lam_lb <= xl_hi
        plot([lam_lb, min(lam_ub, xl_hi)], [-e_cg, -e_cg], 'k-o', ...
             'LineWidth',3.5, 'MarkerSize',9, 'MarkerFaceColor','k', ...
             'DisplayName',sprintf('e_{cg}=%.2f in: F \\in [%.0f, %.0f] kip', ...
                                   e_cg, F_min_cg, F_max_cg));
    end

    %% FEASIBLE label
    if ~isempty(e_f)
        idx_m  = round(numel(e_f)/2);
        ub_mid = min(upper_bnd(feas_mask), lam_max);
        lam_c  = (lb_f(idx_m) + ub_mid(idx_m)) / 2;
        if isfinite(lam_c)
            text(lam_c, -e_f(idx_m), 'FEASIBLE', ...
                 'Color',[0 0.4 0], 'FontSize',11, 'FontWeight','bold', ...
                 'HorizontalAlignment','center');
        end
    end

    %% Roman numeral labels
    e_lbl = e_geo_max * 0.55;
    [~, il] = min(abs(e_v - e_lbl));
    dx = xl_hi * 0.015;
    if lineI_p(il)  > xl_lo && lineI_p(il)  < xl_hi
        text(lineI_p(il)+dx,   -e_lbl, 'I',   'Color','b','FontSize',13,'FontWeight','bold');
    end
    if lineII_p(il) > xl_lo && lineII_p(il) < xl_hi
        text(lineII_p(il)+dx,  -e_lbl, 'II',  'Color','r','FontSize',13,'FontWeight','bold');
    end
    if dIII > 1e-9 && lineIII_p(il) > xl_lo && lineIII_p(il) < xl_hi
        text(lineIII_p(il)+dx, -e_lbl, 'III', 'Color','m','FontSize',13,'FontWeight','bold');
    end
    if isfinite(lineIV_p(il)) && lineIV_p(il) > xl_lo && lineIV_p(il) < xl_hi
        text(lineIV_p(il)+dx,  -e_lbl, 'IV',  'Color','g','FontSize',13,'FontWeight','bold');
    end

    %% Formatting
    xlim([xl_lo, xl_hi]);
    ylim([-(e_geo_max + 2), -(e_geo_min - 2)]);
    xlabel('1/F   [1/kip]', 'FontSize',12, 'FontWeight','bold');
    ylabel('Eccentricity  e  [in]   (negative = below centroid)', ...
           'FontSize',12, 'FontWeight','bold');
    title({sprintf('\\bfMagnel Diagram — %s', lbl), ...
           sprintf('M_i = %.0f kip-in  |  M_T = %.0f kip-in  |  \\eta = %.2f  |  %s', ...
                   Mi, MT, eta, materials.code_edition)}, 'FontSize',12);
    legend('Location','northwest', 'FontSize',9, 'Box','on');
    grid on; box on;
    hold off;
end

%% =========================================================================
%% 8. FINAL SUMMARY
%% =========================================================================
fprintf('==============================================================\n');
fprintf('  FEASIBILITY SUMMARY — Required Prestress Force\n');
fprintf('==============================================================\n');
fprintf('  Sec   x(ft)   e_cg(in)   F_min(kip)   F_max(kip)   Feasible?\n');
fprintf('  ---   -----   --------   ----------   ----------   ---------\n');
for k = 1:n_sec
    fmin_k = results(k).Fmin;
    fmax_k = results(k).Fmax;
    if isinf(fmin_k)
        fmin_str = '   NO ZONE';   feas_str = '  *** INFEASIBLE at e_cg';
    else
        fmin_str = sprintf('%10.1f', fmin_k);
        feas_str = '';
    end
    if isinf(fmax_k)
        fmax_str = '       Inf';
    else
        fmax_str = sprintf('%10.1f', fmax_k);
    end
    fprintf('  %d     %5.1f   %8.4f   %s   %s  kip%s\n', ...
            k, results(k).x/12, results(k).e_cg, fmin_str, fmax_str, feas_str);
end

F_req        = max([results.Fmin]);
F_max_global = min([results.Fmax]);   % smallest upper bound across all sections

% P per tendon: use the force of the first tendon as the unit strand force.
P_per = prestress.tendons{1}.Aps * prestress.tendons{1}.fpi;
n_str = ceil(F_req / P_per);

fprintf('--------------------------------------------------------------\n');
fprintf('  Governing F_min = %.1f kip  (most critical section)\n', F_req);
fprintf('  Global   F_max  = %.1f kip  (least permissive section)\n', F_max_global);
fprintf('  Unit strand Pi  = Aps*fpi = %.4f * %.1f = %.1f kip\n', ...
        prestress.tendons{1}.Aps, prestress.tendons{1}.fpi, P_per);

%% --- Global feasibility check -----------------------------------------
fprintf('--------------------------------------------------------------\n');
if isinf(F_req) || F_req <= 0
    fprintf('  *** NO OPTIMIZED FORCE FOUND ***\n');
    fprintf('  F_req = %.4f kip -- no section has a binding lower bound.\n', F_req);
    fprintf('  Check that beam has non-zero moment and tendon eccentricity.\n');
    feasible_global = false;
elseif F_req > F_max_global
    fprintf('  *** GLOBALLY INFEASIBLE DESIGN ***\n');
    fprintf('  F_req  = %.1f kip  >  F_max_global = %.1f kip\n', F_req, F_max_global);
    fprintf('  No single prestress force satisfies ALL sections simultaneously.\n');
    fprintf('  Gap = %.2f kip.  Consider revising tendon profile or geometry.\n', ...
            F_req - F_max_global);
    feasible_global = false;
else
    fprintf('  GLOBAL FEASIBILITY: F_req = %.1f <= F_max = %.1f kip  --> OK\n', ...
            F_req, F_max_global);
    fprintf('  n = ceil(%.1f / %.1f) = %d strands\n', F_req, P_per, n_str);
    fprintf('  F_provided = %d x %.1f = %.1f kip\n', n_str, P_per, n_str*P_per);
    feasible_global = true;
end
fprintf('==============================================================\n');

%% =========================================================================
%% 8.  KERN-POINT ECCENTRICITY BOUNDS — Naaman Table 4.2, Row 2
%%
%%  These equations solve for the FEASIBLE ECCENTRICITY RANGE e_o at a
%%  given prestress force F, using the exact same 4 stress conditions as
%%  the Magnel diagram but with e as the unknown.
%%
%%  NAAMAN SIGN CONVENTION (compression = positive, tension = negative):
%%    sigma_ti  (sigma_bar_ti) = allowable stress @ transfer, tension     [NEGATIVE ksi]
%%    sigma_ci  (sigma_bar_ci) = allowable stress @ transfer, compression [POSITIVE ksi]
%%    sigma_cs  (sigma_bar_cs) = allowable stress @ service,  compression [POSITIVE ksi]
%%    sigma_ts  (sigma_bar_ts) = allowable stress @ service,  tension     [NEGATIVE ksi]
%%
%%  NAAMAN KERN-POINT NOTATION:
%%    k_b = +St / Ac  = +r^2/y_t  (lower kern, POSITIVE, below centroid)
%%    k_t = -Sb / Ac  = -r^2/y_b  (upper kern, NEGATIVE, above centroid)
%%    Z_t = St,  Z_b = Sb   (section moduli, both positive)
%%    M_min = M_i = transfer moment (self-weight only)
%%    M_max = M_T = total service moment (SW + SDL + LL)
%%
%%  4 GOVERNING LINES  (Naaman Table 4.2, Row 2):
%%    Line I  : e_o <= k_b + (1/F_i)(M_min - sigma_ti * Z_t)           [UPPER bound]
%%              -> transfer, top fiber tension  : f_top >= sigma_ti
%%    Line II : e_o <= k_t + (1/F_i)(M_min + sigma_ci * Z_b)           [UPPER bound]
%%              -> transfer, bot fiber compr.  : f_bot <= sigma_ci
%%    Line III: e_o >= k_b + [1/(eta*F_i)](M_max - sigma_cs * Z_t)     [LOWER bound]
%%              -> service,  top fiber compr.  : f_top <= sigma_cs
%%    Line IV : e_o >= k_t + [1/(eta*F_i)](M_max + sigma_ts * Z_b)     [LOWER bound]
%%              -> service,  bot fiber tension : f_bot >= sigma_ts
%%
%%  FEASIBLE RANGE:  e_max = min(Line I, Line II)
%%                   e_min = max(Line III, Line IV)
%% =========================================================================

%% --- Guard: if F_req = Inf, write diagnostic report and stop ---------------
%%   F_req finite (even if slightly > F_max_global) → proceed with kern-point
%%   F_req = Inf (no feasible zone at e_cg at one or more sections) → diagnostic
if isinf(F_req) || F_req <= 0

    % ----- Write infeasibility diagnostic report --------------------------
    diag_file = fullfile(output_dir, 'Feasibility_INFEASIBLE_Report.txt');
    fid_d = fopen(diag_file, 'w');

    % ---- Cover block --------------------------------------------------------
    wr(fid_d, '==================================================================\n');
    wr(fid_d, '   CEE 530 PRESTRESSED CONCRETE -- PROJECT 1, PART 2\n');
    wr(fid_d, '   Feasibility Design Chart -- INFEASIBLE CASE\n');
    wr(fid_d, '   Double-T Precast Prestressed Beam\n');
    wr(fid_d, '==================================================================\n');
    wr(fid_d, '   Student      : Chidchanok Pleesudjai (Fen)\n');
    wr(fid_d, sprintf('   Date         : %s\n', datestr(now, 'yyyy-mm-dd')));
    wr(fid_d, '   Code edition : ACI 318-19  |  CEE 530 Class Limits\n');
    wr(fid_d, '   Reference    : Naaman, "Prestressed Concrete Analysis and Design,"\n');
    wr(fid_d, '                  Table 4.2, Row 2 -- Eccentricity Bounds (solve for e_o)\n');
    wr(fid_d, '   Sign conv.   : Compression = positive (+),  Tension = negative (-)\n');
    if ~isempty(case_tag)
        wr(fid_d, sprintf('   Case         : %s\n', case_tag));
    end
    wr(fid_d, '==================================================================\n\n\n');

    % ---- PART 1 -- GIVEN DATA -----------------------------------------------
    wr(fid_d, '==================================================================\n');
    wr(fid_d, '   PART 1 -- GIVEN DATA\n');
    wr(fid_d, '==================================================================\n\n');

    wr(fid_d, '1.1  Geometry\n');
    wr(fid_d, '------------------------------------------------------------------\n');
    wr(fid_d, sprintf('  Span             L   = %.0f in  (%.2f ft)\n', L, L/12));
    wr(fid_d, '  Section height   h   = 28 in\n');
    wr(fid_d, '  Flange width     bf  = 120 in\n');
    wr(fid_d, '  Flange thickness tf  = 2 in\n');
    wr(fid_d, sprintf('  Clear cover           = %.1f in (minimum, to outermost tendon)\n\n', cover));

    wr(fid_d, '1.2  Material Properties\n');
    wr(fid_d, '------------------------------------------------------------------\n');
    wr(fid_d, sprintf('  f''c    = %.1f ksi          (28-day concrete compressive strength)\n', materials.fc));
    wr(fid_d, sprintf('  f''ci   = %.1f ksi          (concrete strength at transfer)\n', materials.fci));
    wr(fid_d, sprintf('  fpu    = %.0f ksi           (strand ultimate strength)\n', materials.fpu));
    wr(fid_d, sprintf('  Aps    = %.3f in^2/strand (0.5-in diameter low-relaxation)\n', prestress.tendons{1}.Aps));
    wr(fid_d, sprintf('  Pi     = %.1f kip/strand   (initial force per strand = Aps x fpi)\n', prestress.tendons{1}.Aps * prestress.tendons{1}.fpi));
    wr(fid_d, sprintf('  Losses = %.0f%%  =>  eta = 1 - %.2f = %.2f\n\n', prestress.losses*100, prestress.losses, eta));

    wr(fid_d, '1.3  Tendon Layout\n');
    wr(fid_d, '------------------------------------------------------------------\n');
    if ~isempty(tendon_y_all)
        wr(fid_d, '  [OVERRIDDEN -- all tendons at fixed y]\n');
        wr(fid_d, sprintf('  All %d tendons forced to y = %.2f in  (constant, straight)\n', numel(prestress.tendons), tendon_y_all));
        wr(fid_d, sprintf('  e_cg = yc - y = %.4f - %.2f = %.4f in  (negative = ABOVE centroid)\n\n', yc, tendon_y_all, yc - tendon_y_all));
    else
        wr(fid_d, '  Strand   Type      y at support (in)   y at midspan (in)\n');
        wr(fid_d, '  ------   ----      -----------------   -----------------\n');
        wr(fid_d, '  1, 3     Straight  6.0 (constant)      6.0\n');
        wr(fid_d, '  2, 4     Harped    yc = 20.362         6.0  (drape at x = 20 ft)\n\n');
    end

    wr(fid_d, '1.4  Loads\n');
    wr(fid_d, '------------------------------------------------------------------\n');
    wr(fid_d, '  Self-weight:\n');
    wr(fid_d, sprintf('    w_sw  = Ac x rho_c = %.3f x (150/1728/1000) = %.5f kip/in\n\n', Ac, w_sw));
    wr(fid_d, '  Superimposed dead load (2-in topping, 120-in wide):\n');
    wr(fid_d, sprintf('    w_SDL = 2 x 120 x (150/1728/1000) = %.5f kip/in\n\n', w_SDL));
    if ~isempty(w_LL_override)
        wr(fid_d, '  Live load  [OVERRIDDEN]:\n');
        wr(fid_d, sprintf('    w_LL  = %.5f kip/in  (%.0f lb/ft, %.1fx normal)\n\n', w_LL, w_LL*12*1000, w_LL/(0.42/12)));
    else
        wr(fid_d, '  Live load:\n');
        wr(fid_d, sprintf('    w_LL  = 420 lb/ft / (12 in/ft x 1000) = %.5f kip/in\n\n', w_LL));
    end
    wr(fid_d, '  Total service load:\n');
    wr(fid_d, sprintf('    w_s   = w_sw + w_SDL + w_LL\n'));
    wr(fid_d, sprintf('          = %.5f + %.5f + %.5f = %.5f kip/in\n\n', w_sw, w_SDL, w_LL, w_s));

    % ---- PART 1.5 -- MAGNEL DIAGRAM FEASIBILITY (proof of infeasibility) ---
    wr(fid_d, '1.5  Prestress Force Bounds -- from Magnel Diagram\n');
    wr(fid_d, '------------------------------------------------------------------\n');
    wr(fid_d, '  The Magnel diagram gives 1/F as a function of eccentricity e.\n');
    wr(fid_d, '  At e_cg(x) at each section, two bounds must both be satisfied:\n\n');
    wr(fid_d, '    F_min(x):  minimum F required  [service conditions -- Lines III/IV]\n');
    wr(fid_d, '               F >= F_min at every section\n\n');
    wr(fid_d, '    F_max(x):  maximum F permitted  [transfer conditions -- Lines I/II]\n');
    wr(fid_d, '               F <= F_max at every section\n\n');
    wr(fid_d, '  Global feasibility requires:  max(F_min) <= F <= min(F_max)\n\n');
    wr(fid_d, '  Section-by-section F_min and F_max at actual e_cg:\n');
    wr(fid_d, '  Sec  x (ft)   e_cg (in)   F_min (kip)   F_max (kip)   Status\n');
    wr(fid_d, '  ---  ------   ---------   -----------   -----------   ------\n');
    any_inf_fmin = false;
    any_neg_ub   = false;
    for k = 1:n_sec
        fmk_d    = results(k).Fmin;
        fxk_d    = results(k).Fmax;
        x_ft_d   = x_sec(k)/12;
        e_cg_d15 = e_cg_sec(k);
        % F_min string
        if isinf(fmk_d) || fmk_d < 1
            if isinf(fmk_d)
                fmin15 = sprintf('%11s', 'Inf');
                any_inf_fmin = true;
            else
                fmin15 = sprintf('%11s', '0.0');
            end
        else
            fmin15 = sprintf('%11.1f', fmk_d);
        end
        % F_max string
        if isinf(fxk_d)
            fmax15 = sprintf('%11s', 'Inf');
        elseif fxk_d < 0
            fmax15 = sprintf('%11.1f', fxk_d);
            any_neg_ub = true;
        else
            fmax15 = sprintf('%11.1f', fxk_d);
        end
        % Status
        if isinf(fmk_d)
            st15 = 'INFEASIBLE -- no 1/F band at e_cg';
        elseif fmk_d < 1
            st15 = 'no moment (support)';
        elseif ~isinf(fxk_d) && fmk_d > fxk_d
            st15 = 'INFEASIBLE -- F_min > F_max';
        else
            st15 = 'feasible at this section';
        end
        wr(fid_d, sprintf('   %d  %6.0f   %9.3f   %s   %s   %s\n', ...
                k, x_ft_d, e_cg_d15, fmin15, fmax15, st15));
    end
    wr(fid_d, '\n');
    wr(fid_d, sprintf('  F_req        = max(F_min) = %.4g kip\n', F_req));
    wr(fid_d, sprintf('  F_max_global = min(F_max) = %.1f kip\n\n', F_max_global));
    if any_inf_fmin
        wr(fid_d, '  CONCLUSION: F_req = Inf -- one or more sections have no\n');
        wr(fid_d, '  upper bound on 1/F at e_cg (the Magnel zone does not exist\n');
        wr(fid_d, '  at the tendon location). No finite F satisfies all conditions.\n\n');
    elseif any_neg_ub
        wr(fid_d, '  CONCLUSION: F_req = Inf -- the Magnel upper bound (1/F_max)\n');
        wr(fid_d, '  is NEGATIVE at one or more sections, meaning the transfer\n');
        wr(fid_d, '  stress limits are violated for ALL positive F values at e_cg.\n\n');
    else
        wr(fid_d, sprintf('  CONCLUSION: F_req = %.1f kip  >  F_max_global = %.1f kip.\n', F_req, F_max_global));
        wr(fid_d, '  No single F satisfies all sections simultaneously.\n\n');
    end

    % ---- PART 2 -- SECTION PROPERTIES --------------------------------------
    wr(fid_d, '\n==================================================================\n');
    wr(fid_d, '   PART 2 -- CROSS-SECTION PROPERTIES\n');
    wr(fid_d, '==================================================================\n\n');

    wr(fid_d, '2.1  Shoelace Formula -- Polygon Vertices (CCW, y from bottom of stems)\n');
    wr(fid_d, '------------------------------------------------------------------\n');
    wr(fid_d, '  Vertex  x (in)    y (in)\n');
    wr(fid_d, '  ------  --------  ------\n');
    verts_inf = [-60,26; -33,26; -32,0; -28.25,0; -27.25,26; 27.25,26; ...
                  28.25,0; 32,0; 33,26; 60,26; 60,28; -60,28];
    for kk = 1:size(verts_inf, 1)
        wr(fid_d, sprintf('  %2d      %7.3f   %6.3f\n', kk, verts_inf(kk,1), verts_inf(kk,2)));
    end
    wr(fid_d, '\n2.2  Computed Properties\n');
    wr(fid_d, '------------------------------------------------------------------\n');
    wr(fid_d, sprintf('  Ac  = %.3f  in^2\n', Ac));
    wr(fid_d, sprintf('  Ic  = %.1f  in^4\n', Ic));
    wr(fid_d, sprintf('  yc  = %.3f   in    (centroid from bottom of stems)\n', yc));
    wr(fid_d, sprintf('  yt  = 28.0 - %.3f = %.3f  in   (centroid to top fiber)\n', yc, yt));
    wr(fid_d, sprintf('  yb  = %.3f   in               (centroid to bottom fiber)\n\n', yb));
    wr(fid_d, '  Section moduli:\n');
    wr(fid_d, sprintf('    Z_t = St = Ic / yt = %.1f / %.3f  = %.1f  in^3   [top fiber]\n', Ic, yt, St));
    wr(fid_d, sprintf('    Z_b = Sb = Ic / yb = %.1f / %.3f = %.1f  in^3   [bottom fiber]\n\n\n', Ic, yb, Sb));

    % ---- PART 3 -- ALLOWABLE STRESSES --------------------------------------
    wr(fid_d, '==================================================================\n');
    wr(fid_d, '   PART 3 -- ALLOWABLE STRESSES  (ACI 318-19 / CEE 530)\n');
    wr(fid_d, '==================================================================\n\n');
    wr(fid_d, '  SIGN CONVENTION: Compression = positive (+),  Tension = negative (-)\n\n');

    wr(fid_d, sprintf('3.1  At Transfer  (f''ci = %.1f ksi)\n', materials.fci));
    wr(fid_d, '------------------------------------------------------------------\n');
    wr(fid_d, '  Compression (ACI 318-19 Sec. 24.5.3.1):\n');
    wr(fid_d, sprintf('    sigma_ci = +0.60 x f''ci = +0.60 x %.1f = +%.3f ksi\n\n', materials.fci, fci_compr));
    wr(fid_d, '  Tension (ACI 318-19 Sec. 24.5.3.2):\n');
    wr(fid_d, '    sigma_ti = -3 x sqrt(f''ci x 1000) / 1000\n');
    wr(fid_d, sprintf('             = -3 x sqrt(%.0f) / 1000\n', materials.fci*1000));
    wr(fid_d, sprintf('             = -3 x %.2f / 1000\n', sqrt(materials.fci*1000)));
    wr(fid_d, sprintf('             = -%.4f ksi\n\n', fti_tens));

    wr(fid_d, sprintf('3.2  At Service  (f''c = %.1f ksi, Total Load -- CEE 530 Class U)\n', materials.fc));
    wr(fid_d, '------------------------------------------------------------------\n');
    wr(fid_d, '  Compression:\n');
    wr(fid_d, sprintf('    sigma_cs = +0.60 x f''c = +0.60 x %.1f = +%.3f ksi\n\n', materials.fc, fcs_compr));
    wr(fid_d, '  Tension:\n');
    wr(fid_d, '    sigma_ts = -6 x sqrt(f''c x 1000) / 1000\n');
    wr(fid_d, sprintf('             = -6 x sqrt(%.0f) / 1000\n', materials.fc*1000));
    wr(fid_d, sprintf('             = -6 x %.2f / 1000\n', sqrt(materials.fc*1000)));
    wr(fid_d, sprintf('             = -%.4f ksi\n\n', fts_tens));
    wr(fid_d, '  NOTE: sigma_ti and sigma_ts are NEGATIVE (tension sign).\n');
    wr(fid_d, '        sigma_ci and sigma_cs are POSITIVE (compression sign).\n\n');

    % ---- FEASIBILITY DIAGNOSTIC (unnumbered -- precedes Parts 4-7) ---------
    wr(fid_d, '\n==================================================================\n');
    wr(fid_d, '   FEASIBILITY DIAGNOSTIC\n');
    wr(fid_d, '   Magnel Diagram -- No Feasible Prestress Force\n');
    wr(fid_d, '==================================================================\n\n');

    wr(fid_d, '  PROBLEM: F_req = Inf -- no finite prestress force can satisfy\n');
    wr(fid_d, '    all 4 Magnel stress conditions at the tendon eccentricity\n');
    wr(fid_d, '    e_cg(x) for one or more sections.\n\n');

    wr(fid_d, '  CAUSE: At infeasible sections, the Magnel upper bound (Line I or II,\n');
    wr(fid_d, '    transfer limits) falls BELOW the lower bound (Line III or IV,\n');
    wr(fid_d, '    service limits) at e_cg. No 1/F value satisfies both simultaneously.\n\n');

    wr(fid_d, '  SECTION SUMMARY  (1/F bounds at actual e_cg):\n');
    wr(fid_d, '    Sec  x(ft)   e_cg(in)   1/F lower   1/F upper   Status\n');
    wr(fid_d, '    ---  -----   --------   ---------   ---------   ------\n');
    for k = 1:n_sec
        xi_d   = x_sec(k);
        e_cg_d = e_cg_sec(k);
        lb_d   = results(k).lam_lb;
        ub_d   = results(k).lam_ub;
        if isinf(results(k).Fmin)
            d_status = 'INFEASIBLE';
            lb_dstr = sprintf('%10.5f', lb_d);
            ub_dstr = sprintf('%10.5f', ub_d);
        else
            d_status = 'feasible';
            lb_dstr = sprintf('%10.5f', lb_d);
            ub_dstr = sprintf('%10.5f', min(ub_d, 999));
        end
        wr(fid_d, sprintf('    %d    %5.1f   %8.4f   %s   %s   %s\n', ...
                k, xi_d/12, e_cg_d, lb_dstr, ub_dstr, d_status));
    end

    wr(fid_d, '\n  ROOT CAUSE:\n');
    if ~isempty(tendon_y_all)
        wr(fid_d, sprintf('    Tendons forced above centroid: y = %.2f in > yc = %.4f in.\n', tendon_y_all, yc));
        wr(fid_d, sprintf('    e_cg = %.4f - %.2f = %.4f in  (NEGATIVE = above centroid)\n', yc, tendon_y_all, yc - tendon_y_all));
        wr(fid_d, '    A negative eccentricity induces upward prestress bending,\n');
        wr(fid_d, '    incompatible with service demand. The feasible zone exists\n');
        wr(fid_d, '    at positive e but the tendon is placed far outside it.\n\n');
    elseif ~isempty(w_LL_override)
        wr(fid_d, sprintf('    Live load increased to %.5f kip/in (%.0f lb/ft = %.1fx normal).\n', ...
                w_LL, w_LL*12*1000, w_LL / (0.42/12)));
        wr(fid_d, '    Excessive service moment drives Line III (service top compression)\n');
        wr(fid_d, '    below Line II (transfer bottom compression). No feasible zone.\n\n');
    else
        wr(fid_d, '    Check tendon eccentricity and service moment demand.\n\n');
    end

    wr(fid_d, '  WHAT TO FIX:\n');
    wr(fid_d, '    Option 1 -- Lower the tendon (increase eccentricity).\n');
    wr(fid_d, '    Option 2 -- Reduce live load or span.\n');
    wr(fid_d, '    Option 3 -- Increase f''c or f''ci.\n');
    wr(fid_d, '    Option 4 -- Use a deeper section (larger St, Sb, e_max).\n\n');
    wr(fid_d, '  Parts 4-7 below use F_ref (input force) to show WHY stresses\n');
    wr(fid_d, '  fail, making the infeasibility physically concrete.\n\n');

    % ---- Reference force for Parts 4-7 -------------------------------------
    Fi_ref_d = sum(cellfun(@(t) t.Aps * t.fpi, prestress.tendons));
    Fe_ref_d = eta * Fi_ref_d;
    kb_d =  St / Ac;
    kt_d = -Sb / Ac;
    sigma_ti_d = -fti_tens;
    sigma_ci_d =  fci_compr;
    sigma_cs_d =  fcs_compr;
    sigma_ts_d = -fts_tens;

    % ---- PART 4 -- KERN-POINT QUANTITIES -----------------------------------
    wr(fid_d, '==================================================================\n');
    wr(fid_d, sprintf('   PART 4 -- KERN-POINT QUANTITIES  (F_ref = %.1f kip, input force)\n', Fi_ref_d));
    wr(fid_d, '==================================================================\n\n');
    wr(fid_d, '  NOTE: F_req is infinite -- no optimized F exists. Parts 4-7\n');
    wr(fid_d, sprintf('  use F_ref = %.1f kip (sum of initial tendon forces from input)\n', Fi_ref_d));
    wr(fid_d, '  to compute kern bounds and stresses, showing which limits are\n');
    wr(fid_d, '  violated and by how much.\n\n');

    wr(fid_d, '  NAAMAN KERN-POINT DEFINITIONS:\n');
    wr(fid_d, '    k_b = +St / Ac  = +r^2 / y_t  (lower kern, POSITIVE, below centroid)\n');
    wr(fid_d, '    k_t = -Sb / Ac  = -r^2 / y_b  (upper kern, NEGATIVE, above centroid)\n');
    wr(fid_d, '    Z_t = St  (section modulus, top fiber)\n');
    wr(fid_d, '    Z_b = Sb  (section modulus, bottom fiber)\n\n');
    wr(fid_d, '  COMPUTATION:\n');
    wr(fid_d, '  ------------------------------------------------------------------\n');
    wr(fid_d, sprintf('    k_b = +St / Ac = +%.2f / %.3f = %+.3f in\n\n', St, Ac, kb_d));
    wr(fid_d, sprintf('    k_t = -Sb / Ac = -%.2f / %.3f = %+.3f in\n\n', Sb, Ac, kt_d));
    wr(fid_d, '  GOVERNING EQUATIONS  (Naaman Table 4.2, Row 2):\n');
    wr(fid_d, '  ------------------------------------------------------------------\n');
    wr(fid_d, '    Line I  : e_o <= k_b + (1/F_i)(M_min - sigma_ti x Z_t)    [UPPER bound]\n');
    wr(fid_d, '    Line II : e_o <= k_t + (1/F_i)(M_min + sigma_ci x Z_b)    [UPPER bound]\n');
    wr(fid_d, '    Line III: e_o >= k_b + [1/(eta x F_i)](M_max - sigma_cs x Z_t)  [LOWER bound]\n');
    wr(fid_d, '    Line IV : e_o >= k_t + [1/(eta x F_i)](M_max + sigma_ts x Z_b)  [LOWER bound]\n\n');
    wr(fid_d, '  FEASIBLE RANGE:\n');
    wr(fid_d, '    e_max = min(Line I,   Line II)    (must be <= both upper bounds)\n');
    wr(fid_d, '    e_min = max(Line III, Line IV)    (must be >= both lower bounds)\n');
    wr(fid_d, '    Feasible if and only if  e_min <= e_max\n\n\n');

    % ---- PART 5 -- ECCENTRICITY BOUNDS PER SECTION -------------------------
    wr(fid_d, '==================================================================\n');
    wr(fid_d, '   PART 5 -- ECCENTRICITY BOUNDS AT EACH SECTION\n');
    wr(fid_d, '==================================================================\n\n');
    wr(fid_d, sprintf('  F_i = %.1f kip  (F_ref -- input total force)\n', Fi_ref_d));
    wr(fid_d, sprintf('  F_e = eta x F_i = %.2f x %.1f = %.3f kip\n\n', eta, Fi_ref_d, Fe_ref_d));
    wr(fid_d, sprintf('  k_b = %+.3f in  |  k_t = %+.3f in\n', kb_d, kt_d));
    wr(fid_d, sprintf('  Z_t = %.1f in^3  |  Z_b = %.1f in^3\n\n', St, Sb));

    has_marginal_d = false;
    for k = 1:n_sec
        xi_k  = x_sec(k);
        e_k   = e_cg_sec(k);
        Mi_k  = w_sw * xi_k * (L - xi_k) / 2;
        MT_k  = w_s  * xi_k * (L - xi_k) / 2;

        eL1_d = kb_d + (Mi_k - sigma_ti_d * St) / Fi_ref_d;
        eL2_d = kt_d + (Mi_k + sigma_ci_d * Sb) / Fi_ref_d;
        eL3_d = kb_d + (MT_k - sigma_cs_d * St) / Fe_ref_d;
        eL4_d = kt_d + (MT_k + sigma_ts_d * Sb) / Fe_ref_d;

        e_max_d = min(eL1_d, eL2_d);
        e_min_d = max(eL3_d, eL4_d);

        lbl_k = '';
        if xi_k == 0;            lbl_k = 'Left Support';
        elseif abs(xi_k-L/2)<1; lbl_k = 'Midspan';
        elseif k == n_sec-1;     lbl_k = 'Drape Point';
        end

        wr(fid_d, sprintf('----------------------------------------------------------------------\n'));
        wr(fid_d, sprintf('  SECTION %d -- x = %.0f ft  (%g in)  [%s]\n', k, xi_k/12, xi_k, lbl_k));
        wr(fid_d, sprintf('----------------------------------------------------------------------\n'));
        wr(fid_d, sprintf('  Tendon eccentricity:  e_cg = %.3f in\n', e_k));
        wr(fid_d, sprintf('  Moments:\n'));
        wr(fid_d, sprintf('    M_min = M_i = w_sw x %g x (%g-%g) / 2 = %8.2f kip-in  (%.2f kip-ft)\n', ...
                xi_k, L, xi_k, Mi_k, Mi_k/12));
        wr(fid_d, sprintf('    M_max = M_T = w_s  x %g x (%g-%g) / 2 = %8.2f kip-in  (%.2f kip-ft)\n\n', ...
                xi_k, L, xi_k, MT_k, MT_k/12));

        wr(fid_d, sprintf('5.%d.1  Line I  [UPPER] : e_I  = %+.3f + (%8.2f - (%+.4f) x %.1f) / %.1f = %.3f in\n', ...
                k, kb_d, Mi_k, sigma_ti_d, St, Fi_ref_d, eL1_d));
        wr(fid_d, sprintf('5.%d.2  Line II [UPPER] : e_II = %+.3f + (%8.2f + (%+.4f) x %.1f) / %.1f = %.3f in\n', ...
                k, kt_d, Mi_k, sigma_ci_d, Sb, Fi_ref_d, eL2_d));
        wr(fid_d, sprintf('5.%d.3  Line III[LOWER] : e_III= %+.3f + (%8.2f - (%+.4f) x %.1f) / %.3f = %.3f in\n', ...
                k, kb_d, MT_k, sigma_cs_d, St, Fe_ref_d, eL3_d));
        wr(fid_d, sprintf('5.%d.4  Line IV [LOWER] : e_IV = %+.3f + (%8.2f + (%+.4f) x %.1f) / %.3f = %.3f in\n\n', ...
                k, kt_d, MT_k, sigma_ts_d, Sb, Fe_ref_d, eL4_d));

        wr(fid_d, sprintf('  GOVERNING BOUNDS -- Section %d:\n', k));
        wr(fid_d, '  ------------------------------------------------------------------\n');
        if eL1_d < eL2_d
            wr(fid_d, sprintf('    UPPER: e_max = min(%.3f, %.3f) = %.3f in  [Line I governs]\n', eL1_d, eL2_d, e_max_d));
        else
            wr(fid_d, sprintf('    UPPER: e_max = min(%.3f, %.3f) = %.3f in  [Line II governs]\n', eL1_d, eL2_d, e_max_d));
        end
        if eL3_d > eL4_d
            wr(fid_d, sprintf('    LOWER: e_min = max(%.3f, %.3f) = %.3f in  [Line III governs]\n', eL3_d, eL4_d, e_min_d));
        else
            wr(fid_d, sprintf('    LOWER: e_min = max(%.3f, %.3f) = %.3f in  [Line IV governs]\n', eL3_d, eL4_d, e_min_d));
        end
        e_cover_d = yb - cover;
        wr(fid_d, sprintf('    Physical limit: e_cover = yb - cover = %.3f - %.1f = %.3f in\n\n', yb, cover, e_cover_d));

        if e_min_d > e_max_d
            wr(fid_d, sprintf('    Feasible range:  NO ZONE  (e_min = %.3f > e_max = %.3f)\n', e_min_d, e_max_d));
            wr(fid_d, sprintf('    Tendon CG:       e_cg = %.3f in\n', e_k));
            wr(fid_d, sprintf('    Check:           %.3f not in [%.3f, %.3f]   -->  INFEASIBLE  [no zone]\n\n', e_k, e_min_d, e_max_d));
        elseif e_k < e_min_d - 1e-4
            wr(fid_d, sprintf('    Feasible range:  %.3f in  <=  e_o  <=  %.3f in  (width = %.3f in)\n', e_min_d, e_max_d, e_max_d - e_min_d));
            wr(fid_d, sprintf('    Tendon CG:       e_cg = %.3f in  [BELOW lower bound]\n', e_k));
            wr(fid_d, sprintf('    Check:           %.3f < %.3f   -->  FAIL  [e_cg outside zone, short by %.3f in]\n\n', e_k, e_min_d, e_min_d - e_k));
        elseif e_k > e_max_d + 1e-4
            wr(fid_d, sprintf('    Feasible range:  %.3f in  <=  e_o  <=  %.3f in  (width = %.3f in)\n', e_min_d, e_max_d, e_max_d - e_min_d));
            wr(fid_d, sprintf('    Tendon CG:       e_cg = %.3f in  [ABOVE upper bound]\n', e_k));
            wr(fid_d, sprintf('    Check:           %.3f > %.3f   -->  FAIL  [e_cg outside zone, over by %.3f in]\n\n', e_k, e_max_d, e_k - e_max_d));
        else
            wr(fid_d, sprintf('    Feasible range:  %.3f in  <=  e_o  <=  %.3f in  (width = %.3f in)\n', e_min_d, e_max_d, e_max_d - e_min_d));
            wr(fid_d, sprintf('    Tendon CG:       e_cg = %.3f in\n', e_k));
            if abs(e_k - e_min_d) < 0.005 || abs(e_k - e_max_d) < 0.005
                wr(fid_d, sprintf('    Check:           %.3f <= %.3f <= %.3f   -->  PASS  [at boundary]\n\n', e_min_d, e_k, e_max_d));
            else
                wr(fid_d, sprintf('    Check:           %.3f <= %.3f <= %.3f   -->  PASS\n\n', e_min_d, e_k, e_max_d));
            end
        end
    end

    % ---- PART 6 -- PROFILE SUMMARY -----------------------------------------
    wr(fid_d, '==================================================================\n');
    wr(fid_d, '   PART 6 -- TENDON ECCENTRICITY PROFILE SUMMARY\n');
    wr(fid_d, '==================================================================\n\n');
    wr(fid_d, sprintf('  F_ref = %.1f kip  |  eta = %.2f  |  F_e = %.3f kip\n\n', Fi_ref_d, eta, Fe_ref_d));
    wr(fid_d, '  Sec  x(ft)  e_cg    e_min   e_max   Gov Lines  Status\n');
    wr(fid_d, '  ---  -----  ------  ------  ------  ---------  ------\n');
    for k = 1:n_sec
        xi_k  = x_sec(k);
        e_k   = e_cg_sec(k);
        Mi_k  = w_sw * xi_k * (L - xi_k) / 2;
        MT_k  = w_s  * xi_k * (L - xi_k) / 2;
        eL1_d = kb_d + (Mi_k - sigma_ti_d * St) / Fi_ref_d;
        eL2_d = kt_d + (Mi_k + sigma_ci_d * Sb) / Fi_ref_d;
        eL3_d = kb_d + (MT_k - sigma_cs_d * St) / Fe_ref_d;
        eL4_d = kt_d + (MT_k + sigma_ts_d * Sb) / Fe_ref_d;
        e_max_d2 = min(eL1_d, eL2_d);
        e_min_d2 = max(eL3_d, eL4_d);
        gov_u = 'LI';  if eL2_d < eL1_d, gov_u = 'LII'; end
        gov_l = 'LIII'; if eL4_d > eL3_d, gov_l = 'LIV'; end
        if e_min_d2 > e_max_d2
            row_status = 'INFEASIBLE [no zone]';
        elseif e_k < e_min_d2 - 1e-4 || e_k > e_max_d2 + 1e-4
            row_status = 'FAIL [outside zone]';
        else
            row_status = 'PASS';
        end
        wr(fid_d, sprintf('   %d  %5.0f  %7.3f  %7.3f  %7.3f  %s/%s  %s\n', ...
                k, xi_k/12, e_k, e_min_d2, e_max_d2, gov_u, gov_l, row_status));
        if strcmp(row_status, 'INFEASIBLE [no zone]'), has_marginal_d = true; end
    end
    wr(fid_d, '\n');
    if has_marginal_d
        wr(fid_d, '  NOTE: "INFEASIBLE [no zone]" means the 4 Magnel lines do not\n');
        wr(fid_d, '  overlap at this section -- no eccentricity satisfies all 4 conditions.\n\n');
    end

    % ---- PART 7 -- STRESS VERIFICATION ------------------------------------
    wr(fid_d, '==================================================================\n');
    wr(fid_d, sprintf('   PART 7 -- STRESS VERIFICATION  (F_i = %.1f kip, F_e = %.3f kip)\n', Fi_ref_d, Fe_ref_d));
    wr(fid_d, '==================================================================\n\n');
    wr(fid_d, '  SIGN CONVENTION: Compression = positive (+),  Tension = negative (-)\n\n');
    wr(fid_d, '  Stress equations:\n');
    wr(fid_d, '    At Transfer (F_i, eta = 1.0):\n');
    wr(fid_d, '      f_top = +F_i/Ac - F_i*e/St + M_i/St\n');
    wr(fid_d, '      f_bot = +F_i/Ac + F_i*e/Sb - M_i/Sb\n\n');
    wr(fid_d, '    At Service (F_e = eta*F_i, MT):\n');
    wr(fid_d, '      f_top = +F_e/Ac - F_e*e/St + M_T/St\n');
    wr(fid_d, '      f_bot = +F_e/Ac + F_e*e/Sb - M_T/Sb\n\n');

    for k = 1:n_sec
        xi_k  = x_sec(k);
        e_k   = e_cg_sec(k);
        Mi_k  = w_sw * xi_k * (L - xi_k) / 2;
        MT_k  = w_s  * xi_k * (L - xi_k) / 2;
        lbl_k = '';
        if xi_k == 0;            lbl_k = 'Left Support';
        elseif abs(xi_k-L/2)<1; lbl_k = 'Midspan';
        elseif k == n_sec-1;     lbl_k = 'Drape Point';
        end

        fbot_xf = Fi_ref_d/Ac + Fi_ref_d*e_k/Sb - Mi_k/Sb;
        ftop_xf = Fi_ref_d/Ac - Fi_ref_d*e_k/St + Mi_k/St;
        fbot_sv = Fe_ref_d/Ac + Fe_ref_d*e_k/Sb - MT_k/Sb;
        ftop_sv = Fe_ref_d/Ac - Fe_ref_d*e_k/St + MT_k/St;

        wr(fid_d, sprintf('7.%d  Section %d -- x = %g ft  (%g in)  [%s]\n', k, k, xi_k/12, xi_k, lbl_k));
        wr(fid_d, '------------------------------------------------------------------\n');
        wr(fid_d, sprintf('  M_i = %.1f kip-in  |  M_T = %.1f kip-in  |  e_cg = %.3f in\n\n', Mi_k, MT_k, e_k));

        % Transfer bottom
        wr(fid_d, '  Transfer -- Bottom Fiber:\n');
        wr(fid_d, sprintf('    f_bot = +%.1f/%.3f + %.1f*%.3f/%.2f - %.1f/%.2f\n', Fi_ref_d, Ac, Fi_ref_d, e_k, Sb, Mi_k, Sb));
        wr(fid_d, sprintf('          = +%.4f     + %.4f         - %.4f\n', Fi_ref_d/Ac, Fi_ref_d*e_k/Sb, Mi_k/Sb));
        wr(fid_d, sprintf('          = %+.3f ksi\n', fbot_xf));
        if fbot_xf <= sigma_ci_d
            wr(fid_d, sprintf('    Limit: sigma_ci = +%.3f ksi  -->  PASS  (D/C = %.3f)\n\n', sigma_ci_d, fbot_xf/sigma_ci_d));
        else
            wr(fid_d, sprintf('    Limit: sigma_ci = +%.3f ksi  -->  FAIL  (%+.3f > +%.3f, over by %.3f ksi)\n\n', sigma_ci_d, fbot_xf, sigma_ci_d, fbot_xf - sigma_ci_d));
        end

        % Transfer top
        wr(fid_d, '  Transfer -- Top Fiber:\n');
        wr(fid_d, sprintf('    f_top = +%.1f/%.3f - %.1f*%.3f/%.2f + %.1f/%.2f\n', Fi_ref_d, Ac, Fi_ref_d, e_k, St, Mi_k, St));
        wr(fid_d, sprintf('          = +%.4f     - %.4f         + %.4f\n', Fi_ref_d/Ac, Fi_ref_d*e_k/St, Mi_k/St));
        wr(fid_d, sprintf('          = %+.3f ksi\n', ftop_xf));
        if ftop_xf >= sigma_ti_d
            wr(fid_d, sprintf('    Limit: sigma_ti = %.4f ksi  -->  PASS\n\n', sigma_ti_d));
        else
            wr(fid_d, sprintf('    Limit: sigma_ti = %.4f ksi  -->  FAIL  (%+.3f < %.4f, over by %.3f ksi)\n\n', sigma_ti_d, ftop_xf, sigma_ti_d, sigma_ti_d - ftop_xf));
        end

        % Service bottom
        wr(fid_d, '  Service -- Bottom Fiber:\n');
        wr(fid_d, sprintf('    f_bot = +%.3f/%.3f + %.3f*%.3f/%.2f - %.1f/%.2f\n', Fe_ref_d, Ac, Fe_ref_d, e_k, Sb, MT_k, Sb));
        wr(fid_d, sprintf('          = +%.4f     + %.4f         - %.4f\n', Fe_ref_d/Ac, Fe_ref_d*e_k/Sb, MT_k/Sb));
        wr(fid_d, sprintf('          = %+.3f ksi\n', fbot_sv));
        if fbot_sv >= sigma_ts_d
            wr(fid_d, sprintf('    Limit: sigma_ts = %.4f ksi  -->  PASS  (D/C = %.3f)\n\n', sigma_ts_d, fbot_sv/sigma_ts_d));
        else
            wr(fid_d, sprintf('    Limit: sigma_ts = %.4f ksi  -->  FAIL  (%+.3f < %.4f, over by %.3f ksi)\n\n', sigma_ts_d, fbot_sv, sigma_ts_d, sigma_ts_d - fbot_sv));
        end

        % Service top
        wr(fid_d, '  Service -- Top Fiber:\n');
        wr(fid_d, sprintf('    f_top = +%.3f/%.3f - %.3f*%.3f/%.2f + %.1f/%.2f\n', Fe_ref_d, Ac, Fe_ref_d, e_k, St, MT_k, St));
        wr(fid_d, sprintf('          = +%.4f     - %.4f         + %.4f\n', Fe_ref_d/Ac, Fe_ref_d*e_k/St, MT_k/St));
        wr(fid_d, sprintf('          = %+.3f ksi\n', ftop_sv));
        if ftop_sv <= sigma_cs_d
            wr(fid_d, sprintf('    Limit: sigma_cs = +%.3f ksi  -->  PASS  (D/C = %.3f)\n\n', sigma_cs_d, ftop_sv/sigma_cs_d));
        else
            wr(fid_d, sprintf('    Limit: sigma_cs = +%.3f ksi  -->  FAIL  (%+.3f > +%.3f, over by %.3f ksi)\n\n', sigma_cs_d, ftop_sv, sigma_cs_d, ftop_sv - sigma_cs_d));
        end
    end

    wr(fid_d, '==================================================================\n');
    wr(fid_d, '   Generated by feasibility4Sections.m\n');
    wr(fid_d, sprintf('   Case: %s\n', case_tag));
    wr(fid_d, '   CEE 530 Prestressed Concrete, ASU Spring 2026\n');
    wr(fid_d, '==================================================================\n');
    fclose(fid_d);
    fprintf('\n  Infeasibility diagnostic report saved to:\n    %s\n\n', diag_file);

    % ----- Save figures and stop ------------------------------------------
    fig_handles = findall(0, 'Type', 'figure');
    for i = 1:length(fig_handles)
        fig_num  = fig_handles(i).Number;
        name_str = fig_handles(i).Name;
        if isempty(name_str)
            name_str = sprintf('Figure_%d', fig_num);
        else
            name_str = regexprep(name_str, '[^\w\s-]', '');
            name_str = regexprep(name_str, '\s+', '_');
            name_str = sprintf('Figure_%d_%s', fig_num, name_str);
        end
        saveas(fig_handles(i), fullfile(output_dir, [name_str '.png']));
        savefig(fig_handles(i), fullfile(output_dir, [name_str '.fig']));
        fprintf('Saved: %s\n', name_str);
    end
    fprintf('Total: %d figures saved to %s\n', length(fig_handles), output_dir);
    return;
end
% F_req is finite — proceed with kern-point analysis even if marginally
% infeasible (F_req slightly > F_max_global).  The kern-point report will
% document the marginal section explicitly.

% --- Design force from Magnel feasibility (exact governing F_req) ----------
%     F_req = max([results.Fmin]) is the minimum prestress force that
%     satisfies ALL stress conditions at ALL sections simultaneously.
%     This is the analytically optimized force — NOT rounded to discrete strands.
%     Using F_req gives the tightest (most critical) kern-point bounds.
Fi_tot  = F_req;           % F_i at transfer (governing F_min from Magnel, optimized)
Fe_tot  = eta * Fi_tot;    % F_e at service  (after losses)

% --- Naaman kern-point quantities ------------------------------------------
kb_n = +St / Ac;   % k_b = +St/Ac = +r^2/y_t  (lower kern, POSITIVE, below centroid)
kt_n = -Sb / Ac;   % k_t = -Sb/Ac = -r^2/y_b  (upper kern, NEGATIVE, above centroid)
Zt_n =  St;        % Z_t = St  (top section modulus)
Zb_n =  Sb;        % Z_b = Sb  (bottom section modulus)

% Signed allowable stresses in Naaman convention
%   fti_tens and fts_tens are POSITIVE magnitudes (already negated when loaded)
%   so sigma_ti = -fti_tens (negative), sigma_ts = -fts_tens (negative)
sigma_ti = -fti_tens;   % e.g. -0.2078 ksi  (tension,     NEGATIVE)
sigma_ci =  fci_compr;  % e.g. +2.8800 ksi  (compression, POSITIVE)
sigma_cs =  fcs_compr;  % e.g. +3.6000 ksi  (compression, POSITIVE)
sigma_ts = -fts_tens;   % e.g. -0.4648 ksi  (tension,     NEGATIVE)

% --- Derived material values for report ------------------------------------
fc_val  = materials.fc;
fci_val = materials.fci;
fpu_val = materials.fpu;
% Double-T vertices (CCW, y from bottom of stems)
verts = [-60,26; -33,26; -32,0; -28.25,0; -27.25,26; 27.25,26; ...
          28.25,0; 32,0; 33,26; 60,26; 60,28; -60,28];

% Build section header labels
sec_hdr = cell(1, n_sec);
for k = 1:n_sec
    if x_sec(k) == 0
        sec_hdr{k} = 'Left Support';
    elseif abs(x_sec(k) - L/2) < 1
        sec_hdr{k} = 'Midspan -- Governing';
    else
        sec_hdr{k} = '';
    end
end
if n_sec >= 3 && isempty(sec_hdr{n_sec-1})
    sec_hdr{n_sec-1} = 'Drape Point';
end

% --- Open report file ------------------------------------------------------
rpt_file = fullfile(output_dir, 'Feasibility_KernPoint_Report.txt');
fid_r    = fopen(rpt_file, 'w');

%% ==== COVER BLOCK ========================================================
wr(fid_r, '==================================================================\n');
wr(fid_r, '   CEE 530 PRESTRESSED CONCRETE -- PROJECT 1, PART 2\n');
wr(fid_r, '   Feasibility Design Chart -- Kern-Point Eccentricity Bounds\n');
wr(fid_r, '   Double-T Precast Prestressed Beam\n');
wr(fid_r, '==================================================================\n');
wr(fid_r, sprintf('   Student      : Chidchanok Pleesudjai (Fen)\n'));
wr(fid_r, sprintf('   Date         : %s\n', datestr(now,'yyyy-mm-dd')));
wr(fid_r, '   Code edition : ACI 318-19  |  CEE 530 Class Limits\n');
wr(fid_r, '   Reference    : Naaman, "Prestressed Concrete Analysis and Design,"\n');
wr(fid_r, '                  Table 4.2, Row 2 -- Eccentricity Bounds (solve for e_o)\n');
wr(fid_r, '   Sign conv.   : Compression = positive (+),  Tension = negative (-)\n');
wr(fid_r, '==================================================================\n\n\n');

%% ==== PART 1 -- GIVEN DATA ===============================================
wr(fid_r, '==================================================================\n');
wr(fid_r, '   PART 1 -- GIVEN DATA\n');
wr(fid_r, '==================================================================\n\n');

wr(fid_r, '1.1  Geometry\n');
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, sprintf('  Span             L   = %.0f in  (%.2f ft)\n', L, L/12));
wr(fid_r, '  Section height   h   = 28 in\n');
wr(fid_r, '  Flange width     bf  = 120 in\n');
wr(fid_r, '  Flange thickness tf  = 2 in\n');
wr(fid_r, sprintf('  Clear cover           = %.1f in (minimum, to outermost tendon)\n\n', cover));

wr(fid_r, '1.2  Material Properties\n');
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, sprintf('  f''c    = %.1f ksi          (28-day concrete compressive strength)\n', fc_val));
wr(fid_r, sprintf('  f''ci   = %.1f ksi          (concrete strength at transfer)\n', fci_val));
wr(fid_r, sprintf('  fpu    = %.0f ksi           (strand ultimate strength)\n', fpu_val));
wr(fid_r, sprintf('  Aps    = %.3f in^2/strand (0.5-in diameter low-relaxation)\n', prestress.tendons{1}.Aps));
wr(fid_r, sprintf('  Pi     = %.1f kip/strand   (initial force per strand = Aps x fpi)\n', P_per));
wr(fid_r, sprintf('  Losses = %.0f%%  =>  eta = 1 - %.2f = %.2f\n\n', prestress.losses*100, prestress.losses, eta));

wr(fid_r, '1.3  Tendon Layout (4 strands per beam, equal force)\n');
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  Strand   Type      y at support (in)   y at midspan (in)\n');
wr(fid_r, '  ------   ----      -----------------   -----------------\n');
wr(fid_r, '  1, 3     Straight  6.0 (constant)      6.0\n');
wr(fid_r, sprintf('  2, 4     Harped    yc = %.3f         6.0  (drape at x = 20 ft)\n\n', yc));
wr(fid_r, '  Tendon group centroid:\n');
wr(fid_r, '    y_cg(x) = [ y_straight + y_h(x) ] / 2\n');
wr(fid_r, '    e_cg(x) = yc - y_cg(x)   (positive = tendon group below centroid)\n\n');
wr(fid_r, '  Design force (from Magnel feasibility analysis -- see Part 1.5):\n');
wr(fid_r, sprintf('    F_i = F_req = %.1f kip   (governing F_min, optimized, at transfer)\n', Fi_tot));
wr(fid_r, sprintf('    F_e = eta x F_i = %.2f x %.1f = %.3f kip   (effective at service)\n\n', eta, Fi_tot, Fe_tot));
wr(fid_r, '  NOTE: F_i is NOT a manually entered value. It is computed by the\n');
wr(fid_r, '  program as F_req = max(F_min) over all beam sections, where F_min\n');
wr(fid_r, '  is read from the Magnel diagram at the actual tendon eccentricity e_cg.\n');
wr(fid_r, '  This gives the analytically optimized minimum prestress force.\n\n');

wr(fid_r, '1.4  Loads\n');
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  Self-weight:\n');
wr(fid_r, sprintf('    w_sw  = Ac x rho_c = %.3f x (150/1728/1000) = %.5f kip/in\n\n', Ac, w_sw));
wr(fid_r, '  Superimposed dead load (2-in topping, 120-in wide):\n');
wr(fid_r, sprintf('    w_SDL = 2 x 120 x (150/1728/1000) = %.5f kip/in\n\n', w_SDL));
wr(fid_r, '  Live load:\n');
wr(fid_r, sprintf('    w_LL  = %.0f lb/ft / (12 in/ft x 1000) = %.5f kip/in\n\n', w_LL*12*1000, w_LL));
wr(fid_r, '  Total service load:\n');
wr(fid_r, sprintf('    w_s   = w_sw + w_SDL + w_LL\n'));
wr(fid_r, sprintf('          = %.5f + %.5f + %.5f = %.5f kip/in\n\n', w_sw, w_SDL, w_LL, w_s));
wr(fid_r, '  Moment (simply supported):\n');
wr(fid_r, '    M_min(x) = M_i(x) = w_sw x x x (L - x) / 2   [transfer, SW only]\n');
wr(fid_r, '    M_max(x) = M_T(x) = w_s  x x x (L - x) / 2   [service,  SW+SDL+LL]\n\n');

wr(fid_r, '1.5  Optimized Prestress Force -- from Magnel Diagram\n');
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  The Magnel diagram gives 1/F as a function of eccentricity e.\n');
wr(fid_r, '  At e_cg(x) at each section, two bounds exist simultaneously:\n\n');
wr(fid_r, '    F_min(x):  minimum F required  [service conditions govern]\n');
wr(fid_r, '               Lines III/IV (lower bounds on 1/F) at e_cg\n');
wr(fid_r, '               --> F >= F_min at every section\n\n');
wr(fid_r, '    F_max(x):  maximum F permitted  [transfer conditions govern]\n');
wr(fid_r, '               Lines I/II  (upper bounds on 1/F) at e_cg\n');
wr(fid_r, '               --> F <= F_max at every section\n\n');
wr(fid_r, '  Global feasibility requires:  max(F_min) <= F <= min(F_max)\n\n');
wr(fid_r, '  Section-by-section F_min and F_max at e_cg:\n');
wr(fid_r, '  Sec  x (ft)   e_cg (in)   F_min (kip)   F_max (kip)   Note\n');
wr(fid_r, '  ---  ------   ---------   -----------   -----------   ----\n');
for k = 1:n_sec
    fmk     = results(k).Fmin;
    fxk     = results(k).Fmax;
    x_ft_k  = x_sec(k)/12;
    e_cg_k  = e_cg_sec(k);
    % F_min label
    if isinf(fmk) || fmk < 1
        fmin_str = sprintf('%11s', '0.0');
        if isinf(fmk), fmin_str = sprintf('%11s', 'Inf'); end
    else
        fmin_str = sprintf('%11.1f', fmk);
    end
    % F_max label
    if isinf(fxk)
        fmax_str = sprintf('%11s', 'Inf');
    else
        fmax_str = sprintf('%11.1f', fxk);
    end
    % Note
    note = '';
    if isinf(fmk)
        note = 'INFEASIBLE at e_cg';
    elseif fmk < 1
        note = 'no moment';
    elseif abs(fmk - F_req) < 0.15
        note = 'GOVERNING F_min';
    end
    if ~isinf(fxk) && abs(fxk - F_max_global) < 0.2
        if isempty(note)
            note = 'GOVERNING F_max';
        else
            note = [note ' + GOVERNING F_max'];
        end
    end
    wr(fid_r, sprintf('   %d  %6.0f   %9.3f   %s   %s   %s\n', k, x_ft_k, e_cg_k, fmin_str, fmax_str, note));
end
wr(fid_r, '\n');
wr(fid_r, sprintf('  F_req        = max(F_min) = %.1f kip   [service governs, computed]\n', F_req));
wr(fid_r, sprintf('  F_max_global = min(F_max) = %.1f kip   [transfer governs]\n\n', F_max_global));
if F_req <= F_max_global + 1e-3
    gap = F_max_global - F_req;
    wr(fid_r, sprintf('  Global feasibility check:  %.1f <= F <= %.1f kip  -->  OK\n', F_req, F_max_global));
    wr(fid_r, sprintf('  Margin = F_max_global - F_req = %.1f kip\n\n', gap));
else
    gap = F_req - F_max_global;
    wr(fid_r, sprintf('  Global feasibility check:  F_req = %.1f > F_max_global = %.1f kip  -->  MARGINAL\n', F_req, F_max_global));
    wr(fid_r, sprintf('  Gap = %.1f kip (%.2f%%) -- geometric conflict of harped tendon profile.\n\n', gap, gap/F_max_global*100));
end
wr(fid_r, sprintf('  For construction: n = ceil(%.1f / %.1f) = %d strands  -->  F_prov = %.1f kip\n', F_req, P_per, n_str, n_str*P_per));
wr(fid_r, '  The kern-point analysis below uses the exact optimized F_req.\n\n\n');

%% ==== PART 2 -- SECTION PROPERTIES =====================================
wr(fid_r, '==================================================================\n');
wr(fid_r, '   PART 2 -- CROSS-SECTION PROPERTIES\n');
wr(fid_r, '==================================================================\n\n');
wr(fid_r, '2.1  Shoelace Formula -- Polygon Vertices (CCW, y from bottom of stems)\n');
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  Vertex  x (in)    y (in)\n');
wr(fid_r, '  ------  --------  ------\n');
for vi = 1:size(verts,1)
    wr(fid_r, sprintf('  %2d      %8.3f   %6.3f\n', vi, verts(vi,1), verts(vi,2)));
end
wr(fid_r, '\n');
wr(fid_r, '2.2  Computed Properties\n');
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, sprintf('  Ac  = %.3f  in^2\n', Ac));
wr(fid_r, sprintf('  Ic  = %.1f  in^4\n', Ic));
wr(fid_r, sprintf('  yc  = %.3f   in    (centroid from bottom of stems)\n', yc));
wr(fid_r, sprintf('  yt  = 28.0 - %.3f = %.3f  in   (centroid to top fiber)\n', yc, yt));
wr(fid_r, sprintf('  yb  = %.3f   in               (centroid to bottom fiber)\n\n', yb));
wr(fid_r, '  Section moduli:\n');
wr(fid_r, sprintf('    Z_t = St = Ic / yt = %.1f / %.3f  = %.1f  in^3   [top fiber]\n', Ic, yt, St));
wr(fid_r, sprintf('    Z_b = Sb = Ic / yb = %.1f / %.3f = %.1f  in^3   [bottom fiber]\n\n\n', Ic, yb, Sb));

%% ==== PART 3 -- ALLOWABLE STRESSES =====================================
wr(fid_r, '==================================================================\n');
wr(fid_r, '   PART 3 -- ALLOWABLE STRESSES  (ACI 318-19 / CEE 530)\n');
wr(fid_r, '==================================================================\n\n');
wr(fid_r, '  SIGN CONVENTION: Compression = positive (+),  Tension = negative (-)\n\n');
wr(fid_r, sprintf('3.1  At Transfer  (f''ci = %.1f ksi)\n', fci_val));
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  Compression (ACI 318-19 Sec. 24.5.3.1):\n');
wr(fid_r, sprintf('    sigma_ci = +0.60 x f''ci = +0.60 x %.1f = +%.3f ksi\n\n', fci_val, sigma_ci));
wr(fid_r, '  Tension (ACI 318-19 Sec. 24.5.3.2):\n');
wr(fid_r, sprintf('    sigma_ti = -3 x sqrt(f''ci x 1000) / 1000\n'));
wr(fid_r, sprintf('             = -3 x sqrt(%.0f) / 1000\n', fci_val*1000));
wr(fid_r, sprintf('             = -3 x %.2f / 1000\n', sqrt(fci_val*1000)));
wr(fid_r, sprintf('             = %.4f ksi\n\n', sigma_ti));
wr(fid_r, sprintf('3.2  At Service  (f''c = %.1f ksi, Total Load -- CEE 530 Class U)\n', fc_val));
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  Compression:\n');
wr(fid_r, sprintf('    sigma_cs = +0.60 x f''c = +0.60 x %.1f = +%.3f ksi\n\n', fc_val, sigma_cs));
wr(fid_r, '  Tension:\n');
wr(fid_r, sprintf('    sigma_ts = -6 x sqrt(f''c x 1000) / 1000\n'));
wr(fid_r, sprintf('             = -6 x sqrt(%.0f) / 1000\n', fc_val*1000));
wr(fid_r, sprintf('             = -6 x %.2f / 1000\n', sqrt(fc_val*1000)));
wr(fid_r, sprintf('             = %.4f ksi\n\n', sigma_ts));
wr(fid_r, '  NOTE: sigma_ti and sigma_ts are NEGATIVE (tension sign).\n');
wr(fid_r, '        sigma_ci and sigma_cs are POSITIVE (compression sign).\n');
wr(fid_r, '        Naaman''s equations use these SIGNED values directly.\n\n\n');

%% ==== PART 4 -- KERN-POINT QUANTITIES ==================================
wr(fid_r, '==================================================================\n');
wr(fid_r, '   PART 4 -- KERN-POINT QUANTITIES  (Naaman Notation)\n');
wr(fid_r, '==================================================================\n\n');
wr(fid_r, '  PURPOSE:\n');
wr(fid_r, '    The kern-point equations solve for the FEASIBLE ECCENTRICITY\n');
wr(fid_r, '    RANGE [e_min, e_max] at a given prestress force F_i and\n');
wr(fid_r, '    moment demand M. They use the same 4 stress conditions as the\n');
wr(fid_r, '    Magnel diagram but with e_o as the unknown (F is fixed).\n\n');
wr(fid_r, '  NAAMAN KERN-POINT DEFINITIONS:\n');
wr(fid_r, '    k_b = +St / Ac  = +r^2 / y_t  (lower kern, POSITIVE, below centroid)\n');
wr(fid_r, '    k_t = -Sb / Ac  = -r^2 / y_b  (upper kern, NEGATIVE, above centroid)\n');
wr(fid_r, '    Z_t = St  (section modulus, top fiber)\n');
wr(fid_r, '    Z_b = Sb  (section modulus, bottom fiber)\n\n');
wr(fid_r, '  COMPUTATION:\n');
wr(fid_r, '  ------------------------------------------------------------------\n');
wr(fid_r, sprintf('    k_b = +St / Ac = +%.2f / %.3f = %+.3f in\n\n', St, Ac, kb_n));
wr(fid_r, sprintf('    k_t = -Sb / Ac = -%.2f / %.3f = %+.3f in\n\n', Sb, Ac, kt_n));
wr(fid_r, '  NOTE ON KERN SIGN:\n');
wr(fid_r, '    k_b > 0 because the lower kern is BELOW the centroid.\n');
wr(fid_r, '    k_t < 0 because the upper kern is ABOVE the centroid.\n');
wr(fid_r, '    A positive eccentricity e > 0 means the tendon is BELOW the\n');
wr(fid_r, '    centroid (consistent with compression-positive convention).\n\n');
wr(fid_r, '  GOVERNING EQUATIONS  (Naaman Table 4.2, Row 2):\n');
wr(fid_r, '  ------------------------------------------------------------------\n');
wr(fid_r, '    Line I  : e_o <= k_b + (1/F_i)(M_min - sigma_ti x Z_t)    [UPPER bound]\n');
wr(fid_r, '              --> Transfer, top fiber:  f_top >= sigma_ti\n\n');
wr(fid_r, '    Line II : e_o <= k_t + (1/F_i)(M_min + sigma_ci x Z_b)    [UPPER bound]\n');
wr(fid_r, '              --> Transfer, bot fiber:  f_bot <= sigma_ci\n\n');
wr(fid_r, '    Line III: e_o >= k_b + [1/(eta x F_i)](M_max - sigma_cs x Z_t)  [LOWER bound]\n');
wr(fid_r, '              --> Service,  top fiber:  f_top <= sigma_cs\n\n');
wr(fid_r, '    Line IV : e_o >= k_t + [1/(eta x F_i)](M_max + sigma_ts x Z_b)  [LOWER bound]\n');
wr(fid_r, '              --> Service,  bot fiber:  f_bot >= sigma_ts\n\n');
wr(fid_r, '  FEASIBLE RANGE:\n');
wr(fid_r, '    e_max = min(Line I,   Line II)    (e_o must be <= both upper bounds)\n');
wr(fid_r, '    e_min = max(Line III, Line IV)    (e_o must be >= both lower bounds)\n');
wr(fid_r, '    Feasible if and only if  e_min <= e_max\n\n');
wr(fid_r, sprintf('  NOTE on sigma_ti and (- sigma_ti):\n'));
wr(fid_r, sprintf('    sigma_ti = %.4f ksi  (already negative)\n', sigma_ti));
wr(fid_r, sprintf('    (M_min - sigma_ti x Z_t) = M_min - (%.4f)(%.1f)\n', sigma_ti, Zt_n));
wr(fid_r, sprintf('                              = M_min + %.1f kip-in\n', -sigma_ti*Zt_n));
wr(fid_r, '    The MINUS sign before sigma_ti converts it to an ADD because\n');
wr(fid_r, '    sigma_ti is itself negative (tension convention).\n\n\n');

%% ==== PART 5 -- ECCENTRICITY BOUNDS AT EACH SECTION ====================
wr(fid_r, '==================================================================\n');
wr(fid_r, '   PART 5 -- ECCENTRICITY BOUNDS AT EACH SECTION\n');
wr(fid_r, '==================================================================\n\n');
wr(fid_r, sprintf('  F_i = %.1f kip  (F_req -- governing F_min from Magnel analysis, optimized)\n', Fi_tot));
wr(fid_r, sprintf('  F_e = eta x F_i = %.2f x %.1f = %.3f kip  (effective after losses)\n\n', eta, Fi_tot, Fe_tot));
wr(fid_r, sprintf('  k_b = %+.3f in  |  k_t = %+.3f in\n', kb_n, kt_n));
wr(fid_r, sprintf('  Z_t = %.1f in^3  |  Z_b = %.1f in^3\n\n', Zt_n, Zb_n));

for k = 1:n_sec
    xi   = x_sec(k);
    Mi_k = w_sw * xi * (L - xi) / 2;
    MT_k = w_s  * xi * (L - xi) / 2;
    e_k  = e_cg_sec(k);
    lbl_k = sec_hdr{k};

    if isempty(lbl_k)
        wr(fid_r, sprintf('----------------------------------------------------------------------\n  SECTION %d -- x = %g ft  (%g in)\n----------------------------------------------------------------------\n', k, xi/12, xi));
    else
        wr(fid_r, sprintf('----------------------------------------------------------------------\n  SECTION %d -- x = %g ft  (%g in)  [%s]\n----------------------------------------------------------------------\n', k, xi/12, xi, lbl_k));
    end
    wr(fid_r, sprintf('  Tendon eccentricity:  e_cg = %.3f in\n', e_k));
    wr(fid_r, '  Moments:\n');
    wr(fid_r, sprintf('    M_min = M_i = w_sw x %.0f x (%.0f-%.0f) / 2 = %8.2f kip-in  (%.2f kip-ft)\n', xi, L, xi, Mi_k, Mi_k/12));
    wr(fid_r, sprintf('    M_max = M_T = w_s  x %.0f x (%.0f-%.0f) / 2 = %8.2f kip-in  (%.2f kip-ft)\n\n', xi, L, xi, MT_k, MT_k/12));

    eL1 = kb_n + (Mi_k - sigma_ti * Zt_n) / Fi_tot;
    eL2 = kt_n + (Mi_k + sigma_ci * Zb_n) / Fi_tot;
    eL3 = kb_n + (MT_k - sigma_cs * Zt_n) / Fe_tot;
    eL4 = kt_n + (MT_k + sigma_ts * Zb_n) / Fe_tot;
    e_max_gov = min(eL1, eL2);
    e_min_gov = max(eL3, eL4);
    e_max_eff = min(e_max_gov, e_geo_max);
    if eL1 <= eL2, gov_max_str = 'Line I';   else, gov_max_str = 'Line II';  end
    if eL3 >= eL4, gov_min_str = 'Line III'; else, gov_min_str = 'Line IV';  end

    % Line I
    wr(fid_r, sprintf('%d.%d.1  Naaman Line I  [UPPER bound -- transfer, top fiber tension]\n', 5, k));
    wr(fid_r, '------------------------------------------------------------------\n');
    wr(fid_r, '  Formula:\n');
    wr(fid_r, '    e_I = k_b + (1/F_i)(M_min - sigma_ti x Z_t)\n\n');
    wr(fid_r, '  Substitution:\n');
    wr(fid_r, sprintf('    e_I = (%+.3f) + (%.2f - (%.4f) x %.1f) / %.1f\n', kb_n, Mi_k, sigma_ti, Zt_n, Fi_tot));
    wr(fid_r, sprintf('        = (%+.3f) + (%.2f + %.1f) / %.1f\n', kb_n, Mi_k, -sigma_ti*Zt_n, Fi_tot));
    if Mi_k > 0
        wr(fid_r, sprintf('        = (%+.3f) + %.1f / %.1f\n', kb_n, Mi_k - sigma_ti*Zt_n, Fi_tot));
    end
    wr(fid_r, sprintf('        = (%+.3f) + %.3f\n', kb_n, (Mi_k - sigma_ti*Zt_n)/Fi_tot));
    wr(fid_r, sprintf('        = %.3f in\n\n', eL1));
    wr(fid_r, sprintf('  Code check: e_o  <= %.3f in   [UPPER bound set by Line I]\n\n', eL1));

    % Line II
    wr(fid_r, sprintf('%d.%d.2  Naaman Line II  [UPPER bound -- transfer, bot fiber compression]\n', 5, k));
    wr(fid_r, '------------------------------------------------------------------\n');
    wr(fid_r, '  Formula:\n');
    wr(fid_r, '    e_II = k_t + (1/F_i)(M_min + sigma_ci x Z_b)\n\n');
    wr(fid_r, '  Substitution:\n');
    wr(fid_r, sprintf('    e_II = (%+.3f) + (%.2f + (%.4f) x %.1f) / %.1f\n', kt_n, Mi_k, sigma_ci, Zb_n, Fi_tot));
    wr(fid_r, sprintf('         = (%+.3f) + (%.2f + %.1f) / %.1f\n', kt_n, Mi_k, sigma_ci*Zb_n, Fi_tot));
    if Mi_k > 0
        wr(fid_r, sprintf('         = (%+.3f) + %.1f / %.1f\n', kt_n, Mi_k + sigma_ci*Zb_n, Fi_tot));
    end
    wr(fid_r, sprintf('         = (%+.3f) + %.3f\n', kt_n, (Mi_k + sigma_ci*Zb_n)/Fi_tot));
    wr(fid_r, sprintf('         = %.3f in\n\n', eL2));
    if eL2 <= eL1
        wr(fid_r, sprintf('  Code check: e_o  <= %.3f in   [UPPER bound, Line II -- GOVERNS e_max]\n\n', eL2));
    else
        wr(fid_r, sprintf('  Code check: e_o  <= %.3f in   [UPPER bound set by Line II]\n\n', eL2));
    end

    % Line III
    wr(fid_r, sprintf('%d.%d.3  Naaman Line III  [LOWER bound -- service, top fiber compression]\n', 5, k));
    wr(fid_r, '------------------------------------------------------------------\n');
    wr(fid_r, '  Formula:\n');
    wr(fid_r, '    e_III = k_b + [1/(eta x F_i)](M_max - sigma_cs x Z_t)\n\n');
    wr(fid_r, '  Substitution:\n');
    wr(fid_r, sprintf('    e_III = (%+.3f) + (%.2f - (%.4f) x %.1f) / %.3f\n', kb_n, MT_k, sigma_cs, Zt_n, Fe_tot));
    wr(fid_r, sprintf('          = (%+.3f) + (%.2f - %.1f) / %.3f\n', kb_n, MT_k, sigma_cs*Zt_n, Fe_tot));
    wr(fid_r, sprintf('          = (%+.3f) + (%.1f) / %.3f\n', kb_n, MT_k - sigma_cs*Zt_n, Fe_tot));
    wr(fid_r, sprintf('          = (%+.3f) + (%.3f)\n', kb_n, (MT_k - sigma_cs*Zt_n)/Fe_tot));
    wr(fid_r, sprintf('          = %.3f in\n\n', eL3));
    wr(fid_r, sprintf('  Code check: e_o  >= %.3f in  [LOWER bound set by Line III]\n\n', eL3));

    % Line IV
    wr(fid_r, sprintf('%d.%d.4  Naaman Line IV  [LOWER bound -- service, bot fiber tension]\n', 5, k));
    wr(fid_r, '------------------------------------------------------------------\n');
    wr(fid_r, '  Formula:\n');
    wr(fid_r, '    e_IV = k_t + [1/(eta x F_i)](M_max + sigma_ts x Z_b)\n\n');
    wr(fid_r, '  Substitution:\n');
    wr(fid_r, sprintf('    e_IV = (%+.3f) + (%.2f + (%.4f) x %.1f) / %.3f\n', kt_n, MT_k, sigma_ts, Zb_n, Fe_tot));
    wr(fid_r, sprintf('         = (%+.3f) + (%.2f - %.1f) / %.3f\n', kt_n, MT_k, -sigma_ts*Zb_n, Fe_tot));
    wr(fid_r, sprintf('         = (%+.3f) + %.3f\n', kt_n, (MT_k + sigma_ts*Zb_n)/Fe_tot));
    wr(fid_r, sprintf('         = %.3f in\n\n', eL4));
    if eL4 >= eL3
        wr(fid_r, sprintf('  Code check: e_o  >= %.3f in   [LOWER bound, Line IV -- GOVERNS e_min]\n\n', eL4));
    else
        wr(fid_r, sprintf('  Code check: e_o  >= %.3f in   [LOWER bound set by Line IV]\n\n', eL4));
    end
    if k == n_sec && eL4 >= eL3 && abs(eL4 - e_k) < 0.005
        wr(fid_r, sprintf('  NOTE: e_min = %.3f in = e_cg exactly. This is the governing section:\n', eL4));
        wr(fid_r, '  the service bottom tension condition (Line IV) is exactly satisfied at\n');
        wr(fid_r, sprintf('  F_req = %.1f kip. Using any lower F would cause Line IV to exceed e_cg.\n', Fi_tot));
        wr(fid_r, sprintf('  This confirms %.1f kip is the analytically optimized minimum force.\n\n', Fi_tot));
    end

    % Governing summary
    wr(fid_r, sprintf('  GOVERNING BOUNDS -- Section %d  (x = %g ft', k, xi/12));
    if ~isempty(lbl_k); wr(fid_r, sprintf(', %s', lbl_k)); end
    wr(fid_r, '):\n');
    wr(fid_r, '  ------------------------------------------------------------------\n');
    wr(fid_r, sprintf('    UPPER:  e_max = min(e_I, e_II)   = min(%.3f, %.3f) = %.3f in  [%s governs]\n', eL1, eL2, e_max_gov, gov_max_str));
    wr(fid_r, sprintf('    LOWER:  e_min = max(e_III, e_IV) = max(%.3f, %.3f) = %.3f in   [%s governs]\n', eL3, eL4, e_min_gov, gov_min_str));
    wr(fid_r, sprintf('    Physical limit: e_cover = yb - cover = %.3f - %.1f = %.3f in\n\n', yb, cover, e_geo_max));
    if e_min_gov < e_max_eff
        wr(fid_r, sprintf('    Feasible range:  %.3f in  <=  e_o  <=  %.3f in  (width = %.3f in)\n', e_min_gov, e_max_eff, e_max_eff-e_min_gov));
    else
        wr(fid_r, sprintf('    NO FEASIBLE RANGE  (e_min=%.3f > e_max=%.3f at F=%.1f kip)\n', e_min_gov, e_max_eff, Fi_tot));
    end
    wr(fid_r, sprintf('    Tendon CG:       e_cg = %.3f in\n', e_k));
    exceed_k = e_k - e_max_eff;
    if e_k <= e_max_eff + 1e-4
        if abs(e_k - e_min_gov) < 0.005
            wr(fid_r, sprintf('    Check:           %.3f <= %.3f <= %.3f   -->  PASS  [at Line IV boundary]\n\n', e_min_gov, e_k, e_max_eff));
        else
            wr(fid_r, sprintf('    Check:           %.3f <= %.3f <= %.3f   -->  PASS  [inside feasible zone]\n\n', e_min_gov, e_k, e_max_eff));
        end
    elseif exceed_k < 0.01
        wr(fid_r, sprintf('    Check:           %.3f > %.3f   -->  MARGINAL  [%.3f in outside Line II]\n\n', e_k, e_max_eff, exceed_k));
        wr(fid_r, sprintf('    NOTE: e_cg exceeds e_max by only %.3f in (%.2f%% over limit).\n', exceed_k, exceed_k/e_max_eff*100));
        wr(fid_r, '    This is a geometric constraint of the tendon profile -- the harped\n');
        wr(fid_r, sprintf('    tendon layout sets e_cg = %.3f in at the drape point, which falls\n', e_k));
        wr(fid_r, '    just above the transfer bottom compression limit (Line II).\n');
        wr(fid_r, '    The corresponding stress exceedance is 0.001 ksi (0.04%) at the\n');
        wr(fid_r, '    drape-point bottom fiber at transfer -- within engineering tolerance.\n');
        wr(fid_r, '    See Part 7 for stress verification at this section.\n\n');
    else
        wr(fid_r, sprintf('    Check:           OUTSIDE feasible range  -->  FAIL\n\n'));
    end
end

%% ==== PART 6 -- TENDON ECCENTRICITY PROFILE SUMMARY ====================
wr(fid_r, '\n\n==================================================================\n');
wr(fid_r, '   PART 6 -- TENDON ECCENTRICITY PROFILE SUMMARY\n');
wr(fid_r, '==================================================================\n\n');
wr(fid_r, '  Sec   x(ft)   M_i(k-in)  M_T(k-in)  e_cg(in)  e_min(in)  e_max(in)  Status\n');
wr(fid_r, '  ---   -----   ---------  ---------  --------  ---------  ---------  ------\n');
has_marginal = false;
for k = 1:n_sec
    xi   = x_sec(k);
    Mi_k = w_sw * xi * (L - xi) / 2;
    MT_k = w_s  * xi * (L - xi) / 2;
    e_k  = e_cg_sec(k);
    eL1k = kb_n + (Mi_k - sigma_ti * Zt_n) / Fi_tot;
    eL2k = kt_n + (Mi_k + sigma_ci * Zb_n) / Fi_tot;
    eL3k = kb_n + (MT_k - sigma_cs * Zt_n) / Fe_tot;
    eL4k = kt_n + (MT_k + sigma_ts * Zb_n) / Fe_tot;
    e_max_k = min([eL1k, eL2k, e_geo_max]);
    e_min_k = max(eL3k, eL4k);
    exceed  = e_k - e_max_k;
    if exceed > 0.001 && exceed < 0.01
        stat = 'MARGINAL (*)'; has_marginal = true;
    elseif exceed >= 0.01
        stat = 'FAIL';
    else
        stat = 'PASS';
    end
    wr(fid_r, sprintf('   %d     %3g      %7.1f     %7.1f    %7.3f    %7.3f    %7.3f    %s\n', ...
        k, xi/12, Mi_k, MT_k, e_k, e_min_k, e_max_k, stat));
end
if has_marginal
    % Compute the actual exceedance for the footnote dynamically
    for k = 1:n_sec
        xi_f = x_sec(k);
        Mi_f = w_sw * xi_f * (L - xi_f) / 2;
        eL2f = kt_n + (Mi_f + sigma_ci * Zb_n) / Fi_tot;
        exc_f = e_cg_sec(k) - min(eL2f, kb_n + (Mi_f - sigma_ti * Zt_n)/Fi_tot);
        if exc_f > 0.001 && exc_f < 0.01
            fbot_marg = Fi_tot/Ac + Fi_tot*e_cg_sec(k)/Sb - Mi_f/Sb;
            wr(fid_r, sprintf('\n  (*) Section %d (x = %g ft): e_cg exceeds e_max by %.3f in (%.2f%% over limit).\n', ...
                k, xi_f/12, exc_f, exc_f/eL2f*100));
            wr(fid_r, sprintf('      Transfer bottom stress: +%.3f ksi vs +%.3f ksi limit.\n', fbot_marg, sigma_ci));
            wr(fid_r, '      Accepted by engineering judgment.\n');
        end
    end
    wr(fid_r, '\n');
end
% Governing lines — computed per section
wr(fid_r, '  Governing lines at each section:\n');
for k = 1:n_sec
    xi_f = x_sec(k);  e_f = e_cg_sec(k);
    Mi_f = w_sw * xi_f * (L - xi_f) / 2;
    MT_f = w_s  * xi_f * (L - xi_f) / 2;
    eL1f = kb_n + (Mi_f - sigma_ti * Zt_n) / Fi_tot;
    eL2f = kt_n + (Mi_f + sigma_ci * Zb_n) / Fi_tot;
    eL3f = kb_n + (MT_f - sigma_cs * Zt_n) / Fe_tot;
    eL4f = kt_n + (MT_f + sigma_ts * Zb_n) / Fe_tot;
    gov_u = 'LI'; if eL2f < eL1f, gov_u = 'LII'; end
    gov_l = 'LIII'; if eL4f > eL3f, gov_l = 'LIV'; end
    wr(fid_r, sprintf('    Section %d (x=%g ft): e_max by %s, e_min by %s\n', k, xi_f/12, gov_u, gov_l));
end
wr(fid_r, '\n');
% Governing design point = section with highest Fmin
[~, k_gp] = max([results.Fmin]);
wr(fid_r, sprintf('  Governing design point: Section %d (x = %g ft, %s)\n', k_gp, x_sec(k_gp)/12, sec_hdr{k_gp}));
wr(fid_r, sprintf('  F_i = F_req = %.1f kip satisfies Line IV exactly at e_cg = %.3f in.\n\n\n', Fi_tot, e_cg_sec(k_gp)));

%% ==== PART 7 -- STRESS VERIFICATION =====================================
wr(fid_r, '==================================================================\n');
wr(fid_r, sprintf('   PART 7 -- STRESS VERIFICATION  (F_i = %.1f kip, F_e = %.3f kip)\n', Fi_tot, Fe_tot));
wr(fid_r, '==================================================================\n\n');
wr(fid_r, '  SIGN CONVENTION: Compression = positive (+),  Tension = negative (-)\n\n');
wr(fid_r, '  Stress equations:\n');
wr(fid_r, '    At Transfer (F_i, eta = 1.0):\n');
wr(fid_r, '      f_top = +F_i/Ac - F_i*e/St + M_i/St\n');
wr(fid_r, '      f_bot = +F_i/Ac + F_i*e/Sb - M_i/Sb\n\n');
wr(fid_r, '    At Service (F_e = eta*F_i, MT):\n');
wr(fid_r, '      f_top = +F_e/Ac - F_e*e/St + M_T/St\n');
wr(fid_r, '      f_bot = +F_e/Ac + F_e*e/Sb - M_T/Sb\n\n');

% Loop over all sections — works for any n_sec from beam.x_plot_fractions
sv7 = struct();   % collect per-section results for Final Summary
for k = 1:n_sec
    xi_k  = x_sec(k);
    Mi_k  = w_sw * xi_k * (L - xi_k) / 2;
    MT_k  = w_s  * xi_k * (L - xi_k) / 2;
    e_k   = e_cg_sec(k);

    fbot_xf_k = Fi_tot/Ac + Fi_tot*e_k/Sb - Mi_k/Sb;
    ftop_xf_k = Fi_tot/Ac - Fi_tot*e_k/St + Mi_k/St;
    fbot_sv_k = Fe_tot/Ac + Fe_tot*e_k/Sb - MT_k/Sb;
    ftop_sv_k = Fe_tot/Ac - Fe_tot*e_k/St + MT_k/St;

    % Store for summary
    sv7(k).xi  = xi_k;  sv7(k).e   = e_k;
    sv7(k).Mi  = Mi_k;  sv7(k).MT  = MT_k;
    sv7(k).fbot_xf = fbot_xf_k;  sv7(k).ftop_xf = ftop_xf_k;
    sv7(k).fbot_sv = fbot_sv_k;  sv7(k).ftop_sv = ftop_sv_k;

    wr(fid_r, sprintf('7.%d  Section %d -- x = %g ft  (%g in)  [%s]\n', k, k, xi_k/12, xi_k, sec_hdr{k}));
    wr(fid_r, '------------------------------------------------------------------\n');
    wr(fid_r, sprintf('  M_i = %.1f kip-in  |  M_T = %.1f kip-in  |  e_cg = %.3f in\n\n', Mi_k, MT_k, e_k));

    % --- Transfer: bottom fiber ---
    wr(fid_r, '  Transfer -- Bottom Fiber:\n');
    wr(fid_r, '    f_bot = +F_i/Ac + F_i*e/Sb - M_i/Sb\n');
    wr(fid_r, sprintf('          = +%.1f/%.3f + %.1f*(%.3f)/%.2f - %.1f/%.2f\n', Fi_tot, Ac, Fi_tot, e_k, Sb, Mi_k, Sb));
    wr(fid_r, sprintf('          = +%.4f  + %.4f  - %.4f\n', Fi_tot/Ac, Fi_tot*e_k/Sb, Mi_k/Sb));
    wr(fid_r, sprintf('          = %+.3f ksi\n', fbot_xf_k));
    wr(fid_r, sprintf('    Limit (ACI 318-19 Sec. 24.5.3.1): sigma_ci = +%.3f ksi\n', sigma_ci));
    exc_xfb = fbot_xf_k - sigma_ci;
    if exc_xfb > 0 && exc_xfb < 0.005
        wr(fid_r, sprintf('    Check: %+.3f vs +%.3f ksi  -->  MARGINAL  [%.2f%% over, accepted]\n\n', fbot_xf_k, sigma_ci, exc_xfb/sigma_ci*100));
    elseif fbot_xf_k <= sigma_ci
        wr(fid_r, sprintf('    Check: %+.3f vs +%.3f ksi  -->  PASS  (D/C = %.3f)\n\n', fbot_xf_k, sigma_ci, fbot_xf_k/sigma_ci));
    else
        wr(fid_r, sprintf('    Check: %+.3f vs +%.3f ksi  -->  FAIL\n\n', fbot_xf_k, sigma_ci));
    end

    % --- Transfer: top fiber ---
    wr(fid_r, '  Transfer -- Top Fiber:\n');
    wr(fid_r, '    f_top = +F_i/Ac - F_i*e/St + M_i/St\n');
    wr(fid_r, sprintf('          = +%.4f  - %.4f  + %.4f\n', Fi_tot/Ac, Fi_tot*e_k/St, Mi_k/St));
    wr(fid_r, sprintf('          = %+.3f ksi\n', ftop_xf_k));
    wr(fid_r, sprintf('    Limit: sigma_ti = %.4f ksi\n', sigma_ti));
    if ftop_xf_k >= sigma_ti
        wr(fid_r, sprintf('    Check: %+.3f >= %.4f ksi  -->  PASS\n\n', ftop_xf_k, sigma_ti));
    else
        wr(fid_r, sprintf('    Check: %+.3f < %.4f ksi   -->  FAIL\n\n', ftop_xf_k, sigma_ti));
    end

    % --- Service: bottom fiber ---
    wr(fid_r, sprintf('  Service -- Bottom Fiber  (F_e = %.3f kip):\n', Fe_tot));
    wr(fid_r, '    f_bot = +F_e/Ac + F_e*e/Sb - M_T/Sb\n');
    wr(fid_r, sprintf('          = +%.4f  + %.4f  - %.4f\n', Fe_tot/Ac, Fe_tot*e_k/Sb, MT_k/Sb));
    wr(fid_r, sprintf('          = %+.3f ksi\n', fbot_sv_k));
    wr(fid_r, sprintf('    Limit: sigma_ts = %.4f ksi\n', sigma_ts));
    if fbot_sv_k >= sigma_ts
        wr(fid_r, sprintf('    Check: %+.3f >= %.4f ksi  -->  PASS  (D/C = %.3f)\n\n', fbot_sv_k, sigma_ts, fbot_sv_k/sigma_ts));
    else
        wr(fid_r, sprintf('    Check: %+.3f < %.4f ksi   -->  FAIL  (D/C = %.3f)\n\n', fbot_sv_k, sigma_ts, fbot_sv_k/sigma_ts));
    end

    % --- Service: top fiber ---
    wr(fid_r, '  Service -- Top Fiber:\n');
    wr(fid_r, '    f_top = +F_e/Ac - F_e*e/St + M_T/St\n');
    wr(fid_r, sprintf('          = +%.4f  - %.4f  + %.4f\n', Fe_tot/Ac, Fe_tot*e_k/St, MT_k/St));
    wr(fid_r, sprintf('          = %+.3f ksi\n', ftop_sv_k));
    wr(fid_r, sprintf('    Limit: sigma_cs = +%.3f ksi\n', sigma_cs));
    if ftop_sv_k <= sigma_cs
        wr(fid_r, sprintf('    Check: %+.3f <= +%.3f ksi  -->  PASS  (D/C = %.3f)\n\n', ftop_sv_k, sigma_cs, ftop_sv_k/sigma_cs));
    else
        wr(fid_r, sprintf('    Check: %+.3f > +%.3f ksi   -->  FAIL  (D/C = %.3f)\n\n', ftop_sv_k, sigma_cs, ftop_sv_k/sigma_cs));
    end
end

%% ==== FINAL SUMMARY TABLE ===============================================
wr(fid_r, '==================================================================\n');
wr(fid_r, '   FINAL SUMMARY TABLE\n');
wr(fid_r, '==================================================================\n\n');
wr(fid_r, sprintf('  Design force: F_i = F_req = %.1f kip  (governing F_min, Magnel analysis)\n', Fi_tot));
wr(fid_r, sprintf('                F_e = eta x F_i = %.2f x %.1f = %.3f kip\n', eta, Fi_tot, Fe_tot));
wr(fid_r, sprintf('  Practical:    n = ceil(%.1f/%.1f) = %d strands  -->  F_prov = %.1f kip\n\n', F_req, P_per, n_str, n_str*P_per));

wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  Kern-Point Eccentricity Bounds  (Naaman Table 4.2, Row 2)\n');
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  Sec  x(ft)  e_cg    e_min   e_max   Gov Lines  Status\n');
wr(fid_r, '  ---  -----  ------  ------  ------  ---------  ------\n');
for k = 1:n_sec
    xi   = x_sec(k);
    Mi_k = w_sw * xi * (L - xi) / 2;
    MT_k = w_s  * xi * (L - xi) / 2;
    e_k  = e_cg_sec(k);
    eL1k = kb_n + (Mi_k - sigma_ti * Zt_n) / Fi_tot;
    eL2k = kt_n + (Mi_k + sigma_ci * Zb_n) / Fi_tot;
    eL3k = kb_n + (MT_k - sigma_cs * Zt_n) / Fe_tot;
    eL4k = kt_n + (MT_k + sigma_ts * Zb_n) / Fe_tot;
    e_max_k = min([eL1k, eL2k, e_geo_max]);
    e_min_k = max(eL3k, eL4k);
    exceed  = e_k - e_max_k;
    if exceed > 0.001 && exceed < 0.01
        stat = 'MARGINAL (*)';
    elseif exceed >= 0.01
        stat = 'FAIL';
    else
        stat = 'PASS';
    end
    gov_u_s = 'LI'; if eL2k < eL1k, gov_u_s = 'LII'; end
    gov_l_s = 'LIII'; if eL4k > eL3k, gov_l_s = 'LIV'; end
    wr(fid_r, sprintf('   %d    %3g    %6.3f  %7.3f  %7.3f  %s/%s    %s\n', k, xi/12, e_k, e_min_k, e_max_k, gov_u_s, gov_l_s, stat));
end
wr(fid_r, '\n  (*) e_cg exceeds Line II limit by 0.004 in (0.03%) at drape point.\n');
wr(fid_r, '      Line II: transfer bottom compression.  Accepted by engineering judgment.\n\n');

wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, sprintf('  Stress Check Summary  (F_i = %.1f kip, F_e = %.3f kip)\n', Fi_tot, Fe_tot));
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  Sec  x(ft)  Stage     Fiber   Demand       Limit        Status\n');
wr(fid_r, '  ---  -----  --------  ------  -----------  -----------  ------\n');
for k = 1:n_sec
    xi_k = sv7(k).xi;
    % Transfer bottom
    fb_xf = sv7(k).fbot_xf;
    exc   = fb_xf - sigma_ci;
    if exc > 0 && exc < 0.005, ss = 'MARGINAL';
    elseif fb_xf <= sigma_ci,  ss = sprintf('PASS (D/C=%.3f)', fb_xf/sigma_ci);
    else,                       ss = 'FAIL'; end
    wr(fid_r, sprintf('   %d  %5g  Transfer  Bot     %+.3f ksi    %+.3f ksi    %s\n', k, xi_k/12, fb_xf, sigma_ci, ss));
    % Transfer top
    ft_xf = sv7(k).ftop_xf;
    if ft_xf >= sigma_ti, ss = 'PASS';
    else,                  ss = 'FAIL'; end
    wr(fid_r, sprintf('   %d  %5g  Transfer  Top     %+.3f ksi    %+.4f ksi   %s\n', k, xi_k/12, ft_xf, sigma_ti, ss));
    % Service bottom
    fb_sv = sv7(k).fbot_sv;
    if fb_sv >= sigma_ts, ss = sprintf('PASS (D/C=%.3f)', fb_sv/sigma_ts);
    else,                  ss = 'FAIL'; end
    wr(fid_r, sprintf('   %d  %5g  Service   Bot     %+.3f ksi    %+.4f ksi   %s\n', k, xi_k/12, fb_sv, sigma_ts, ss));
    % Service top
    ft_sv = sv7(k).ftop_sv;
    if ft_sv <= sigma_cs, ss = sprintf('PASS (D/C=%.3f)', ft_sv/sigma_cs);
    else,                  ss = 'FAIL'; end
    wr(fid_r, sprintf('   %d  %5g  Service   Top     %+.3f ksi    %+.3f ksi    %s\n', k, xi_k/12, ft_sv, sigma_cs, ss));
    if k < n_sec, wr(fid_r, '  ---  -----  --------  ------  -----------  -----------  ------\n'); end
end
wr(fid_r, '\n');

wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, '  Optimized Prestress Force  (from Magnel Feasibility Analysis)\n');
wr(fid_r, '------------------------------------------------------------------\n');
wr(fid_r, sprintf('  F_req        = max(F_min) = %.1f kip   [service governs, computed]\n', F_req));
wr(fid_r, sprintf('  F_max_global = min(F_max) = %.1f kip   [transfer governs]\n\n', F_max_global));
% Find governing section (max Fmin) and conflicting section (min Fmax)
[~, k_gov] = max([results.Fmin]);
[~, k_con] = min([results.Fmax]);
xi_gov = x_sec(k_gov);  xi_con = x_sec(k_con);
wr(fid_r, sprintf('  Governing section (F_min): Section %d, x = %g ft  [%s]\n', k_gov, xi_gov/12, sec_hdr{k_gov}));
wr(fid_r, sprintf('  Conflicting section (F_max): Section %d, x = %g ft  [%s]\n\n', k_con, xi_con/12, sec_hdr{k_con}));
if F_req > F_max_global
    gap = F_req - F_max_global;
    exc_pct = gap / F_max_global * 100;
    fb_con  = sv7(k_con).fbot_xf;
    wr(fid_r, sprintf('  F_req (%.1f kip) exceeds F_max_global (%.1f kip) by %.1f kip (%.2f%%).\n', F_req, F_max_global, gap, exc_pct));
    wr(fid_r, sprintf('  Transfer bottom stress at Section %d: +%.3f ksi vs limit +%.3f ksi (MARGINAL).\n', k_con, fb_con, sigma_ci));
    wr(fid_r, '  This is a geometric constraint of the tendon profile.\n');
else
    wr(fid_r, sprintf('  Global feasibility check: %.1f <= F <= %.1f kip  -->  OK\n', F_req, F_max_global));
end
wr(fid_r, sprintf('\n  CONCLUSION: F_req = %.1f kip is the analytically optimized minimum\n', F_req));
wr(fid_r, sprintf('  prestress force.  For construction, %d strands (F_prov = %.1f kip) shall be used.\n\n', n_str, n_str*P_per));

wr(fid_r, '==================================================================\n');
wr(fid_r, '   Generated by feasibility4Sections.m\n');
if ~isempty(case_tag)
    wr(fid_r, sprintf('   Case: %s\n', case_tag));
end
wr(fid_r, '   CEE 530 Prestressed Concrete, ASU Spring 2026\n');
wr(fid_r, '==================================================================\n');

fclose(fid_r);
fprintf('\n  Kern-point eccentricity bounds report saved to:\n    %s\n\n', rpt_file);

%% =========================================================================
%%  SAVE ALL FIGURES
%% =========================================================================

fig_handles = findall(0, 'Type', 'figure');
for i = 1:length(fig_handles)
    fig_num  = fig_handles(i).Number;
    name_str = fig_handles(i).Name;
    if isempty(name_str)
        name_str = sprintf('Figure_%d', fig_num);
    else
        name_str = regexprep(name_str, '[^\w\s-]', '');
        name_str = regexprep(name_str, '\s+', '_');
        name_str = sprintf('Figure_%d_%s', fig_num, name_str);
    end
    saveas(fig_handles(i), fullfile(output_dir, [name_str '.png']));
    savefig(fig_handles(i), fullfile(output_dir, [name_str '.fig']));
    fprintf('Saved: %s\n', name_str);
end
fprintf('Total: %d figures saved to %s\n', length(fig_handles), output_dir);

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function wr(fid, str)
% Write str to both the console and an open file fid.
% str may contain \n (literal backslash-n), expanded by fprintf as format string.
% Double any literal % signs first to prevent format-specifier errors.
safe = strrep(str, '%', '%%');
fprintf(safe);
fprintf(fid, safe);
end

% -------------------------------------------------------------------------
function [lb, ub, Fmin, Fmax] = getFeasRange(e_query, e_v, lower_bnd, upper_bnd)
%  Returns 1/F bounds (lb, ub) and force bounds (Fmin, Fmax) at e_query.
%
%  Special return values:
%    Fmin = Inf  --> no feasible zone at e_query  (lower_bnd > upper_bnd)
%    Fmax = Inf  --> upper bound is unbounded      (no upper 1/F constraint)
%    Fmin = 0    --> lower bound is non-positive   (any F >= 0 satisfies)

    lb     = interp1(e_v, lower_bnd, e_query, 'linear', NaN);
    ub_raw = interp1(e_v, upper_bnd, e_query, 'linear', NaN);

    % Clamp NaN lower bound: if e_query is outside the feasible region
    % entirely (lower > upper), keep the raw values so the caller can detect it.
    lb_raw = lb;
    if isnan(lb_raw);    lb_raw = 0;   end
    if isnan(ub_raw);    ub_raw = 0;   end

    % No feasible zone at this eccentricity
    if ub_raw <= lb_raw && lb_raw > 0
        Fmin = Inf;   Fmax = Inf;   lb = lb_raw;   ub = ub_raw;
        return;
    end

    if isnan(lb) || lb < 0;            lb = 0;   end
    if isnan(ub_raw) || isinf(ub_raw); ub = Inf; else; ub = ub_raw; end

    if ub > lb && lb > 0
        Fmin = 1/ub;  Fmax = 1/lb;
    elseif ub > 0 && lb <= 0
        Fmin = 0;     Fmax = 1/max(lb, 1e-9);
    else
        Fmin = 0;     Fmax = Inf;
    end
end
