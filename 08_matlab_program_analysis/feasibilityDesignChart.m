%% =========================================================================
%%  feasibilityDesignChart.m
%%  Magnel Diagram (Feasibility Design Chart) for Prestressed Beam
%%  Project 1 — Double-T Section, Full Service Load (SW + 2" SDL + LL)
%%
%%  ALGORITHM (Nawy / Lin & Burns notation, compression = POSITIVE):
%%
%%  Fiber stresses (compression +, tension -):
%%    At Transfer (F, M_i = self-weight moment):
%%      f_top_i = F/Ac  -  F·e/St  +  M_i/St
%%      f_bot_i = F/Ac  +  F·e/Sb  -  M_i/Sb
%%    At Service (η·F, M_T = SW + SDL moment):
%%      f_top_s = η·F/Ac  -  η·F·e/St  +  M_T/St
%%      f_bot_s = η·F/Ac  +  η·F·e/Sb  -  M_T/Sb
%%
%%  Allowable stress limits (ACI 318):
%%    Transfer : compression ≤ +fci = 0.60·f'ci
%%               tension    ≥ -fti = -3·√f'ci (ksi)
%%    Service  : compression ≤ +fcs = 0.45·f'c
%%               tension    ≥ -fts = -12·√f'c (ksi)
%%
%%  Magnel lines — Naaman convention (4 governing + 2 additional):
%%    Line I   (top@transfer, TENS):   1/F ≥ (e/St  - 1/Ac) / (fti + M_i/St)   [lower, e>kt]
%%    Line II  (bot@transfer, COMPR):  1/F ≥ (1/Ac  + e/Sb) / (fci + M_i/Sb)  [lower]
%%    Line III (bot@service,  TENS):   1/F ≤ η·(1/Ac + e/Sb) / (M_T/Sb - fts) [upper]
%%    Line IV  (top@service,  COMPR):  1/F ≥ η·(1/Ac - e/St) / (fcs - M_T/St) [lower]
%%    Line V   (top@transfer, COMPR):  1/F ≥ (1/Ac  - e/St) / (fci - M_i/St)  [lower, e<kt, addl]
%%    Line VI  (bot@transfer, TENS):   1/F ≤ (1/Ac  + e/Sb) / (M_i/Sb - fti)  [upper, addl]
%%
%%  Feasible region: max(I,II,IV,V) ≤ 1/F ≤ min(III,VI), with e ≤ e_geo_max
%%
%%  Figure 1 — Magnel Diagram at midspan (e vs 1/F)  ← Naaman axes: e on y, 1/F on x
%%  Figure 2 — Feasibility Zone along beam (e_min/e_max vs x/L)
%%              with tendon CG profile superimposed
%%  Figures auto-saved to:  08_matlab_program_analysis/Project1/
%% =========================================================================
clear; clc; close all;

%% Output folder (all figures saved here as .png and .fig)
out_dir = fullfile(fileparts(mfilename('fullpath')), 'Project1');
if ~exist(out_dir, 'dir');  mkdir(out_dir);  end

fprintf('=================================================================\n');
fprintf(' Prestressed Concrete — Magnel / Feasibility Design Chart\n');
fprintf(' Project 1 | Double-T Section | Full Service Load (SW + SDL + LL)\n');
fprintf('=================================================================\n\n');

%% =========================================================================
%% 1.  CROSS-SECTION GEOMETRY
%%     Origin: y = 0 at BOTTOM of stems, y = 28 at TOP of flange
%%     Vertices listed counter-clockwise (required by shoelace formula)
%% =========================================================================
vertices = [
    -60,    26;   %  1  Bottom-left  of top flange
    -33,    26;   %  2  Inner-top    of left stem (flange/stem junction)
    -32,     0;   %  3  Inner-bottom of left stem
    -28.25,  0;   %  4  Outer-bottom of left stem
    -27.25, 26;   %  5  Outer-top    of left stem
     27.25, 26;   %  6  Outer-top    of right stem
     28.25,  0;   %  7  Outer-bottom of right stem
     32,     0;   %  8  Inner-bottom of right stem
     33,    26;   %  9  Inner-top    of right stem
     60,    26;   % 10  Bottom-right of top flange
     60,    28;   % 11  Top-right    of top flange
    -60,    28;   % 12  Top-left     of top flange
];

[Ac, Ic, yc, yt, yb] = computeSectionProps(vertices);

St = Ic / yt;          % Section modulus, top  fiber [in^3]
Sb = Ic / yb;          % Section modulus, bottom fiber [in^3]
r2 = Ic / Ac;          % Radius of gyration squared [in^2]

fprintf('--- Section Properties ---\n');
fprintf('  Ac  = %9.3f  in^2\n', Ac);
fprintf('  Ic  = %9.1f  in^4\n', Ic);
fprintf('  yc  = %9.4f  in  (centroid from bottom)\n', yc);
fprintf('  yt  = %9.4f  in  (centroid → top  fiber)\n', yt);
fprintf('  yb  = %9.4f  in  (centroid → bot  fiber)\n', yb);
fprintf('  St  = %9.2f  in^3\n', St);
fprintf('  Sb  = %9.2f  in^3\n', Sb);
fprintf('  r²  = %9.4f  in^2\n\n', r2);

%% =========================================================================
%% 2.  MATERIAL PROPERTIES & ALLOWABLE STRESSES
%% =========================================================================
fc   = 6.0;   % 28-day compressive strength f'c  [ksi]
fci  = 4.8;   % Compressive strength at transfer f'ci  [ksi]

%  Compression   (positive)
fci_allow = 0.60 * fci;                        %  2.88  ksi  (ACI 318 §24.5.3.1)
fcs_allow = 0.60 * fc;                         %  3.60  ksi  (CEE530: total load limit 0.60*f'c)

%  Tension  magnitude  (positive, applied as lower limit = −f_t_allow)
fti_allow = 3.0  * sqrt(fci * 1000) / 1000;   %  0.208 ksi  (ACI 318 §24.5.3.2)
fts_allow = 6.0  * sqrt(fc  * 1000) / 1000;   %  0.465 ksi  (CEE530 Class U: 6√f'c)

fprintf('--- Allowable Stresses ---\n');
fprintf('  At Transfer  (f''ci = %.1f ksi):\n', fci);
fprintf('    Compression  fci = +%.4f ksi  (0.60 x f''ci)\n', fci_allow);
fprintf('    Tension      fti = -%.4f ksi  (3*sqrt(f''ci))\n',      fti_allow);
fprintf('  At Service   (f''c  = %.1f ksi):\n', fc);
fprintf('    Compression  fcs = +%.4f ksi  (0.60 x f''c  — CEE530 total load)\n',  fcs_allow);
fprintf('    Tension      fts = -%.4f ksi  (6*sqrt(f''c)  — CEE530 Class U)\n\n',    fts_allow);

%% =========================================================================
%% 3.  BEAM GEOMETRY, LOADS & MOMENTS
%% =========================================================================
L      = (20 + 24 + 20) * 12;       % Simply-supported span  [in] = 768 in
eta    = 1 - 0.15;                   % Prestress efficiency after 15% losses
rho_c  = 150 / 1728 / 1000;          % Concrete density [kip/in^3]  (150 pcf = 8.681e-5 kip/in^3)

w_sw   = Ac * rho_c;                 % Self-weight (uniform) [kip/in]
w_SDL  = 2 * 120 * rho_c;           % 2" concrete topping, 120" wide [kip/in]
w_LL   = 420 / 12 / 1000;           % Live load 420 lb/ft -> kip/in
w_D    = w_sw + w_SDL;               % Total dead load [kip/in]
w_serv = w_sw + w_SDL + w_LL;        % Total service load (SW + SDL + LL) [kip/in]

%  Simply-supported beam:  M(x) = w*x*(L-x)/2
npts   = 1001;
x_vec  = linspace(0, L, npts);       % position along span [in]

M_i_x  = w_sw   .* x_vec .* (L - x_vec) / 2;   % transfer moment  (SW only)       [kip-in]
M_T_x  = w_serv .* x_vec .* (L - x_vec) / 2;   % service  moment  (SW+SDL+LL)     [kip-in]

M_i_mid = w_sw   * L^2 / 8;          % midspan transfer moment  [kip-in]
M_T_mid = w_serv * L^2 / 8;          % midspan service  moment  [kip-in]

fprintf('--- Applied Loads & Moments ---\n');
fprintf('  Span L       = %6.0f in  (%5.1f ft)\n',    L, L/12);
fprintf('  w_sw         = %.6f kip/in  (%.3f kip/ft)\n',  w_sw,  w_sw*12);
fprintf('  w_SDL        = %.6f kip/in  (%.3f kip/ft)\n',  w_SDL, w_SDL*12);
fprintf('  w_LL         = %.6f kip/in  (%.3f kip/ft)\n',  w_LL,  w_LL*12);
fprintf('  w_serv       = %.6f kip/in  (%.3f kip/ft)  [SW+SDL+LL]\n',  w_serv, w_serv*12);
fprintf('  M_i @midspan = %8.2f kip-in  (%7.2f kip-ft)\n', M_i_mid, M_i_mid/12);
fprintf('  M_T @midspan = %8.2f kip-in  (%7.2f kip-ft)\n', M_T_mid, M_T_mid/12);
fprintf('  eta          = %.2f\n\n', eta);

%% =========================================================================
%% 4.  DESIGN POINT  (4 tendons × 25 kip each from input file)
%%     Tendons 1,3: straight,    y = 6 in (from bottom)
%%     Tendons 2,4: trapezoidal, y = 20.36 in at ends → 6 in at midspan
%%     At midspan: all four tendons at y = 6 in → e = yc − 6
%% =========================================================================
F_des  = 4 * 25;                     % Total initial prestress force [kip]
y_t_mid = 6;                         % Tendon CG at midspan [in from bottom]
e_des   = yc - y_t_mid;             % Eccentricity at midspan [in]
lam_des = 1 / F_des;                 % Magnel parameter 1/F [1/kip]

fprintf('--- Design Point (4 × 0.5" strands, Pe = 25 kip each) ---\n');
fprintf('  F_des    = %7.2f  kip\n',   F_des);
fprintf('  e_des    = %7.4f  in   (tendons at y = %.2f in from bottom)\n', e_des, y_t_mid);
fprintf('  1/F_des  = %7.6f  1/kip\n\n', lam_des);

%% =========================================================================
%% 5.  MAGNEL DIAGRAM  (midspan — critical section)
%%     X-axis: eccentricity e [in]
%%     Y-axis: λ = 1/F      [1/kip]
%% =========================================================================
Mi   = M_i_mid;        % scalar, midspan transfer moment
MT   = M_T_mid;        % scalar, midspan service  moment
cover = 2.0;           % clear cover [in]

e_geo_max = yb - cover;          % physical upper eccentricity bound
e_geo_min = -(yt - cover);       % physical lower eccentricity bound (tendon at top)

e_v = linspace(e_geo_min - 2, e_geo_max + 3, 1200);   % eccentricity vector

% ---- Magnel line denominators (Naaman governing convention) ----
% Line I  : Top  @ Transfer — TENSION    limit  [lower bound, binding for e > kt = r²/yt]
dI   = fti_allow + Mi / St;      % > 0 always

% Line II : Bot  @ Transfer — COMPRESSION limit  [lower bound]
dII  = fci_allow + Mi / Sb;      % > 0 always

% Line III: Bot  @ Service  — TENSION    limit  [upper bound]
dIII = MT / Sb   - fts_allow;    % > 0 (service moment creates tension at bottom)

% Line IV : Top  @ Service  — COMPRESSION limit  [lower bound]
dIV  = fcs_allow - MT / St;      % > 0 (service compression within limit)

% Line V  : Top  @ Transfer — COMPRESSION limit  [lower bound, binding for e < kt, additional]
dV   = fci_allow - Mi / St;      % usually > 0

% Line VI : Bot  @ Transfer — TENSION    limit  [upper bound, additional]
dVI  = Mi / Sb   - fti_allow;    % > 0 near midspan (moment large enough)

% ---- Six Magnel lines ----
% Naaman governing (I–IV):
lineI   =       (e_v/St - 1/Ac) ./ dI;    % lower bound (top tension @ transfer, e > kt)
lineII  =       (1/Ac + e_v/Sb) ./ dII;   % lower bound (bot compression @ transfer)
lineIII = eta * (1/Ac + e_v/Sb) ./ dIII;  % upper bound (bot tension @ service)
lineIV  = eta * (1/Ac - e_v/St) ./ dIV;   % lower bound (top compression @ service)

% Additional (V–VI):
if dV > 1e-12
    lineV = (1/Ac - e_v/St) ./ dV;        % lower bound (top compression @ transfer)
else
    lineV = -inf(size(e_v));
end
if dVI > 1e-12
    lineVI = (1/Ac + e_v/Sb) ./ dVI;      % upper bound (bot tension @ transfer)
else
    lineVI = inf(size(e_v));
end

% ---- Envelope ----
lower_bnd = max([lineI; lineII; lineIV; lineV], [], 1);
upper_bnd = min([lineIII; lineVI],              [], 1);

% ---- Feasible zone (λ > 0, upper > lower) ----
feas_mask = (upper_bnd > lower_bnd) & (lower_bnd > 0) & ...
            (e_v <= e_geo_max) & (e_v >= e_geo_min);

% ---- Kern limits ----
k_b = r2 / yb;    % lower kern (max e for no tension at top under centric P)
k_t = r2 / yt;    % upper kern

fprintf('--- Magnel Diagram Parameters (midspan, Naaman convention) ---\n');
fprintf('  dI   (fti+Mi/St) = %+.4f  ksi  [Line I:  top tens @ transfer]\n',  dI);
fprintf('  dII  (fci+Mi/Sb) = %+.4f  ksi  [Line II: bot compr @ transfer]\n', dII);
fprintf('  dIII (MT/Sb-fts) = %+.4f  ksi  [Line III: bot tens @ service]\n',  dIII);
fprintf('  dIV  (fcs-MT/St) = %+.4f  ksi  [Line IV: top compr @ service]\n',  dIV);
fprintf('  Lower kern k_b   = %+.4f  in\n',  k_b);
fprintf('  Upper kern k_t   = %+.4f  in\n\n', k_t);

%% ---- FIGURE 1: Magnel Diagram ----
fig1 = figure('Name','Figure 1 - Magnel Diagram', 'Position',[40 80 920 680]);
hold on;

%  Negate e for plot: class convention e > 0 = above centroid, e < 0 = below centroid
%  (internal math keeps e_old = yc - y_tendon, positive below; we just flip for display)
e_f  = e_v(feas_mask);
lb_f = lower_bnd(feas_mask);
ub_f = upper_bnd(feas_mask);
if ~isempty(e_f)
    fill([lb_f, fliplr(ub_f)], [-e_f, fliplr(-e_f)], ...
         [0.65 0.95 0.65], 'EdgeColor','none', 'FaceAlpha',0.55, ...
         'HandleVisibility','off');
    plot(lb_f, -e_f, 'k-', 'LineWidth',0.5, 'HandleVisibility','off');
    plot(ub_f, -e_f, 'k-', 'LineWidth',0.5, 'HandleVisibility','off');
end

%  Magnel lines — x = 1/F [1/kip],  y = e [in]  (+ above centroid, − below)
hI   = plot(lineI,   -e_v, 'b-',  'LineWidth',2.2, ...
       'DisplayName','I: Top@Transfer (tens. \geq -f_{ti})  [Naaman]');
hII  = plot(lineII,  -e_v, 'r-',  'LineWidth',2.2, ...
       'DisplayName','II: Bot@Transfer (compr. \leq f_{ci})  [Naaman]');
hIII = plot(lineIII, -e_v, 'm--', 'LineWidth',2.2, ...
       'DisplayName','III: Bot@Service  (tens. \geq -f_{ts})');
hIV  = plot(lineIV,  -e_v, 'g-',  'LineWidth',2.2, ...
       'DisplayName','IV: Top@Service  (compr. \leq f_{cs})');
hV   = plot(lineV,   -e_v, 'b:',  'LineWidth',1.4, ...
       'DisplayName','V: Top@Transfer (compr. \leq f_{ci})  [addl]');
hVI  = plot(lineVI,  -e_v, 'r:',  'LineWidth',1.4, ...
       'DisplayName','VI: Bot@Transfer (tens. \geq -f_{ti})  [addl]');

%  Reference lines
xline(0, 'k--', 'LineWidth',0.8, 'HandleVisibility','off');
yline(0, 'k--', 'LineWidth',0.8, 'HandleVisibility','off');
yline(-e_geo_max, 'k:',  'LineWidth',1.8, ...
      'DisplayName',sprintf('e_{cover} = %.2f in  (physical lower bound)', -e_geo_max));
yline(-k_b, 'c-.', 'LineWidth',1.2, ...
      'DisplayName',sprintf('−k_b = %.2f in  (lower kern)', -k_b));
yline( k_t, 'c:',  'LineWidth',1.2, ...
      'DisplayName',sprintf('+k_t = %.2f in  (upper kern)',  k_t));

%  Design point — only plot if it fits within the axis range (otherwise annotate)
%  We will set xlim focused on feasible zone; add arrow annotation for out-of-range point
plot(NaN, NaN, 'rp', 'MarkerSize',14, 'MarkerFaceColor','r', 'LineWidth',1.5, ...
     'DisplayName',sprintf('Design pt (F=%.0f kip) OUTSIDE  -->  1/F=%.5f', F_des, lam_des));

%  Feasible zone label
if ~isempty(e_f)
    idx_mid = round(numel(e_f)/2);
    lam_mid = (lb_f(idx_mid) + ub_f(idx_mid)) / 2;
    text(lam_mid, -e_f(idx_mid), 'FEASIBLE', ...
         'Color',[0.0 0.45 0.0], 'FontSize',12, 'FontWeight','bold', ...
         'HorizontalAlignment','center');
end

%  Line labels — placed at a fixed (negative) e_label
e_label = e_geo_max * 0.6;           % internal e value (positive = below centroid)
[~, il] = min(abs(e_v - e_label));
dx = 0.0003;
text(lineI(il)  + dx, -e_label, 'I',   'Color','b', 'FontSize',12,'FontWeight','bold');
text(lineII(il) + dx, -e_label, 'II',  'Color','r', 'FontSize',12,'FontWeight','bold');
text(lineIII(il)+ dx, -e_label, 'III', 'Color','m', 'FontSize',12,'FontWeight','bold');
text(lineIV(il) + dx, -e_label, 'IV',  'Color','g', 'FontSize',12,'FontWeight','bold');

%  Axes & formatting
xlabel('Magnel parameter,  1/F   [1/kip]', 'FontSize',13,'FontWeight','bold');
ylabel('Eccentricity,  e  [in]   (+ above centroid,  − below centroid)', ...
       'FontSize',13,'FontWeight','bold');
title({'\bfMagnel Diagram  —  Feasibility Design Chart', ...
       sprintf('Double-T Section  |  L = %.0f in (%.1f ft)  |  Full Service Load (SW+SDL+LL)', ...
               L, L/12)}, 'FontSize',13);
legend('Location','northwest','FontSize',9,'Box','on');
grid on;  box on;
ax1 = gca;  ax1.FontSize = 11;

ylim([-(e_geo_max+3), -(e_geo_min-2)]);
if ~isempty(e_f)
    % Focus x-axis on feasible zone — do NOT stretch for out-of-range design point
    xl_lo = max(min(lb_f)*0.6, -0.0002);
    xl_hi = max(ub_f) * 2.5;
    xlim([xl_lo, xl_hi]);
end

% Annotate out-of-range design point with arrow from right edge
ax1_lim = xlim;
ax1_ylim = ylim;
y_star = -e_des;   % -14.36 in
x_arrow_tip  = ax1_lim(2) * 0.98;   % near right edge
x_arrow_tail = ax1_lim(2) * 0.70;
annotation('textarrow', ...
    [x_arrow_tail/ax1_lim(2), x_arrow_tip/ax1_lim(2)], ...
    [(y_star - ax1_ylim(1))/(ax1_ylim(2)-ax1_ylim(1)), ...
     (y_star - ax1_ylim(1))/(ax1_ylim(2)-ax1_ylim(1))], ...
    'String', sprintf('F=%g kip OUTSIDE\\rightarrow (1/F=%.4f)', F_des, lam_des), ...
    'FontSize', 9, 'Color', 'r', 'HeadStyle', 'vback2', ...
    'HeadLength', 8, 'HeadWidth', 6);

hold off;

%  --- Print Magnel line equations for verification ---
fprintf('--- Magnel Line Equations (1/F = a + b*e) — Naaman convention ---\n');
fprintf('  Line I   (top tens @ xfer):  1/F = %+.6g + (%+.6g)*e\n', -1/(Ac*dI),  1/(St*dI));
fprintf('  Line II  (bot compr @ xfer): 1/F = %+.6g + (%+.6g)*e\n',  1/(Ac*dII), 1/(Sb*dII));
fprintf('  Line III (bot tens @ svc):   1/F = %+.6g + (%+.6g)*e\n',  eta/(Ac*dIII), eta/(Sb*dIII));
fprintf('  Line IV  (top compr @ svc):  1/F = %+.6g + (%+.6g)*e\n',  eta/(Ac*dIV), -eta/(St*dIV));
fprintf('\n');

%% =========================================================================
%% 6.  FEASIBILITY ZONE ALONG BEAM  (e_min / e_max vs x/L, for F = F_des)
%%     Solved by rearranging each stress condition for e given F and M(x)
%% =========================================================================
F = F_des;

%  Eccentricity bounds from each of the 8 conditions (vectors over x_vec):
%
%  Transfer (F, M_i_x):
%    Cond I   (top compr ≤ fci):   e ≥ St/Ac + (M_i − fci·St)/F
%    Cond II  (top tens  ≥ −fti):  e ≤ St/Ac + (M_i + fti·St)/F
%    Cond III (bot tens  ≥ −fti):  e ≥ (M_i − fti·Sb)/F − Sb/Ac
%    Cond IV  (bot compr ≤ fci):   e ≤ (M_i + fci·Sb)/F − Sb/Ac
%
%  Service (η·F, M_T_x):
%    Cond V   (top compr ≤ fcs):   e ≥ St/Ac + (M_T − fcs·St)/(η·F)
%    Cond VI  (top tens  ≥ −fts):  e ≤ St/Ac + (M_T + fts·St)/(η·F)
%    Cond VII (bot tens  ≥ −fts):  e ≥ (M_T − fts·Sb)/(η·F) − Sb/Ac
%    Cond VIII(bot compr ≤ fcs):   e ≤ (M_T + fcs·Sb)/(η·F) − Sb/Ac

eF = eta * F;   % effective prestress force at service

e_lb_I    =  St/Ac + (M_i_x - fci_allow*St) / F;       % lower (transfer, top compr)
e_ub_II   =  St/Ac + (M_i_x + fti_allow*St) / F;       % upper (transfer, top tens)
e_lb_III  = (M_i_x - fti_allow*Sb) / F  - Sb/Ac;       % lower (transfer, bot tens)
e_ub_IV   = (M_i_x + fci_allow*Sb) / F  - Sb/Ac;       % upper (transfer, bot compr)

e_lb_V    =  St/Ac + (M_T_x - fcs_allow*St) / eF;      % lower (service, top compr)
e_ub_VI   =  St/Ac + (M_T_x + fts_allow*St) / eF;      % upper (service, top tens)
e_lb_VII  = (M_T_x - fts_allow*Sb) / eF - Sb/Ac;       % lower (service, bot tens) ← governing
e_ub_VIII = (M_T_x + fcs_allow*Sb) / eF - Sb/Ac;       % upper (service, bot compr)

%  Combined envelope
e_min_env = max([e_lb_I;  e_lb_III;  e_lb_V;  e_lb_VII],  [], 1);
e_max_env = min([e_ub_II; e_ub_IV;   e_ub_VI; e_ub_VIII], [], 1);
e_max_env = min(e_max_env, e_geo_max * ones(size(e_max_env)));   % physical limit

%% ---- Tendon CG Profile ----
%  Tendon 1 & 3 (straight, y = 6 in from bottom everywhere)
y_t13 = 6 * ones(size(x_vec));   % [in from bottom]

%  Tendon 2 & 4 (trapezoidal harped):
%    Drape points: x = L/2 ± 144 in  → [0, 240, 528, 768] in
%    y-profile   : [20.36, 6, 6, 20.36] in
xp_h  = [0,  L/2-144,  L/2+144,  L];
yp_h  = [28-7.64,   6,   6,   28-7.64];      %  [20.36, 6, 6, 20.36]
y_t24 = interp1(xp_h, yp_h, x_vec, 'linear');

%  CG of all 4 tendons (2 straight + 2 harped, equal area each)
y_CG  = (y_t13 + y_t24) / 2;    % average = (6 + y_harped)/2
e_CG  = yc - y_CG;              % eccentricity of CG [in]

%  Individual tendon eccentricities
e_t13 = yc - y_t13;             % straight tendons
e_t24 = yc - y_t24;             % harped tendons

%% ---- FIGURE 2: Feasibility Zone ----
fig2 = figure('Name','Figure 2 - Feasibility Zone', 'Position',[980 80 900 620]);
hold on;

x_pct = x_vec / L * 100;   % position as % of span

%  Negate e for display: + above centroid, − below centroid
valid = e_max_env > e_min_env;
fill([x_pct(valid), fliplr(x_pct(valid))], ...
     [-e_min_env(valid), fliplr(-e_max_env(valid))], ...
     [0.65 0.95 0.65], 'EdgeColor','none', 'FaceAlpha',0.50, ...
     'HandleVisibility','off');

%  Eccentricity boundary lines (negated)
%  -e_min_env = upper (tendon closest to centroid allowed by prestress requirement)
%  -e_max_env = lower (cover limit, deepest the tendon can go)
plot(x_pct, -e_min_env, 'b-', 'LineWidth',2.2, ...
     'DisplayName','e  upper bound  (min prestress required)');
plot(x_pct, -e_max_env, 'r-', 'LineWidth',2.2, ...
     'DisplayName','e  lower bound  (cover limit)');

%  Governing condition labels (at midspan location)
[~, ixm] = max(x_pct(x_pct <= 50));   % 50% span index
% Lower-bound conditions: addl(xfer top compr), Naaman Line IV, addl(xfer bot tens), Naaman Line III
gov_min = {'addl(xfer-top-compr)', 'Naaman Line IV', 'addl(xfer-bot-tens)', 'Naaman Line III'};
% Upper-bound conditions: Naaman Line I, Naaman Line II, addl(svc top tens), addl(svc bot compr)
gov_max = {'Naaman Line I', 'Naaman Line II', 'addl(svc-top-tens)', 'addl(svc-bot-compr)'};
lbs = [e_lb_I(ixm); e_lb_V(ixm); e_lb_III(ixm); e_lb_VII(ixm)];
ubs = [e_ub_II(ixm); e_ub_IV(ixm); e_ub_VI(ixm); e_ub_VIII(ixm)];
[~, gi_min] = max(lbs);
[~, gi_max] = min(ubs);
text(52, -e_min_env(ixm)+0.8, ...
     sprintf('governs: %s', gov_min{gi_min}), ...
     'Color','b','FontSize',9,'FontAngle','italic');
text(52, -e_max_env(ixm)-0.8, ...
     sprintf('governs: %s', gov_max{gi_max}), ...
     'Color','r','FontSize',9,'FontAngle','italic');

%  Physical geometric limit
yline(-e_geo_max, 'k:', 'LineWidth',1.8, ...
      'DisplayName',sprintf('e_{cover} = %.2f in  (cover = %.0f in from bottom)', ...
                             -e_geo_max, cover));

%  Individual tendon eccentricities (negated: below centroid → negative)
plot(x_pct, -e_t13, 'b--', 'LineWidth',1.4, ...
     'DisplayName','Tendons 1 & 3 (straight)');
plot(x_pct, -e_t24, 'm--', 'LineWidth',1.4, ...
     'DisplayName','Tendons 2 & 4 (harped)');
plot(x_pct, -e_CG,  'k-',  'LineWidth',2.6, ...
     'DisplayName','Tendon group CG,  e(x)');

%  Reference
yline(0, 'k--', 'LineWidth',0.8, 'HandleVisibility','off');
xline(0,  'k-', 'LineWidth',0.8, 'HandleVisibility','off');
xline(100,'k-', 'LineWidth',0.8, 'HandleVisibility','off');
xline(50,  'k-.','LineWidth',0.6,'HandleVisibility','off');
text(50.5, min(-e_max_env)-1.5, 'Midspan','FontSize',9,'Color',[0.4 0.4 0.4]);

xlabel('Position along span,  x/L  [%]', 'FontSize',13,'FontWeight','bold');
ylabel('Eccentricity,  e  [in]   (+ above centroid,  − below centroid)', ...
       'FontSize',13,'FontWeight','bold');
title({'\bfFeasibility Zone  —  Eccentricity Bounds vs Position', ...
       sprintf('F = %.0f kip  |  \\eta = %.2f  |  L = %.0f in (%.1f ft)  |  Full Service (SW+SDL+LL)', ...
               F, eta, L, L/12)}, 'FontSize',13);
legend('Location','north','FontSize',9,'NumColumns',2,'Box','on');
grid on; box on;
ax2 = gca; ax2.FontSize = 11;
xlim([0, 100]);
ylim([min(-e_max_env)-3, max(-e_min_env)+3]);

hold off;

%% =========================================================================
%% 7.  STRESS VERIFICATION AT MIDSPAN  (design point check)
%% =========================================================================
fprintf('=================================================================\n');
fprintf(' Stress Verification at Midspan\n');
fprintf(' F = %.0f kip | e = %.4f in | η = %.2f\n', F_des, e_des, eta);
fprintf(' M_i = %.2f kip-in | M_T = %.2f kip-in\n', Mi, MT);
fprintf('=================================================================\n');

F = F_des;  e = e_des;

f_top_i = F/Ac  - F*e/St    + Mi/St;
f_bot_i = F/Ac  + F*e/Sb    - Mi/Sb;
f_top_s = eta*F/Ac - eta*F*e/St + MT/St;
f_bot_s = eta*F/Ac + eta*F*e/Sb - MT/Sb;

fprintf('\n  AT TRANSFER (F = %.0f kip):\n', F);
fprintf('    f_top = %+.4f ksi | allow [%+.4f, %+.4f]  →  %s\n', ...
        f_top_i, -fti_allow, fci_allow, chkS(f_top_i,-fti_allow,fci_allow));
fprintf('    f_bot = %+.4f ksi | allow [%+.4f, %+.4f]  →  %s\n', ...
        f_bot_i, -fti_allow, fci_allow, chkS(f_bot_i,-fti_allow,fci_allow));

fprintf('\n  AT SERVICE (η·F = %.0f kip):\n', eta*F);
fprintf('    f_top = %+.4f ksi | allow [%+.4f, %+.4f]  →  %s\n', ...
        f_top_s, -fts_allow, fcs_allow, chkS(f_top_s,-fts_allow,fcs_allow));
fprintf('    f_bot = %+.4f ksi | allow [%+.4f, %+.4f]  →  %s\n', ...
        f_bot_s, -fts_allow, fcs_allow, chkS(f_bot_s,-fts_allow,fcs_allow));

fprintf('\n  Check — design point in Magnel feasible zone:\n');
lam_lo = interp1(e_v, lower_bnd, e_des, 'linear');
lam_hi = interp1(e_v, upper_bnd, e_des, 'linear');
fprintf('    λ_min(e_des) = %.6f 1/kip\n', lam_lo);
fprintf('    λ_des        = %.6f 1/kip\n', lam_des);
fprintf('    λ_max(e_des) = %.6f 1/kip\n', lam_hi);
if lam_des >= lam_lo && lam_des <= lam_hi
    fprintf('    → Design point is WITHIN the feasible zone.  ✓\n');
else
    fprintf('    → Design point is OUTSIDE the feasible zone.  ✗\n');
end

%% =========================================================================
%% 7.5  ECCENTRICITY BOUNDS — KERN-POINT LINEAR EQUATIONS
%%
%%  Naaman kern-point form  (solved for e_o, not 1/F):
%%
%%    k_t = St / Ac     (kern point, top-fiber side,    = r²/yt)
%%    k_b = Sb / Ac     (kern point, bottom-fiber side, = r²/yb)
%%    Z_t = St,  Z_b = Sb   (section moduli)
%%
%%  8 conditions (Naaman textbook labels shown for 4 governing lines):
%%  TRANSFER  (F_i):
%%   [addl]     e ≥  +k_t + (M_i − f_ci·Z_t) / F_i    top compr ≤ +fci  [lower, additional]
%%   [Line I]   e ≤  +k_t + (M_i + f_ti·Z_t) / F_i    top tens  ≥ −fti  [upper bound on e]
%%   [addl]     e ≥  −k_b + (M_i − f_ti·Z_b) / F_i    bot tens  ≥ −fti  [lower, additional]
%%   [Line II]  e ≤  −k_b + (M_i + f_ci·Z_b) / F_i    bot compr ≤ +fci  [upper bound on e]
%%  SERVICE  (η·F):
%%   [Line IV]  e ≥  +k_t + (M_T − f_cs·Z_t) / (η·F)  top compr ≤ +fcs  [lower bound on e]
%%   [addl]     e ≤  +k_t + (M_T + f_ts·Z_t) / (η·F)  top tens  ≥ −fts  [upper, additional]
%%   [Line III] e ≥  −k_b + (M_T − f_ts·Z_b) / (η·F)  bot tens  ≥ −fts  [lower bound on e]
%%   [addl]     e ≤  −k_b + (M_T + f_cs·Z_b) / (η·F)  bot compr ≤ +fcs  [upper, additional]
%%
%%  Naaman 4 governing lines:
%%    e_min = max(Line IV, Line III)   →  lower bound on e
%%    e_max = min(Line I,  Line II)    →  upper bound on e
%%  Additional checks tighten the feasible range further when active.
%% =========================================================================

kt_k = St / Ac;          % kern point: top fiber side  [in]  (= r²/yt)
kb_k = Sb / Ac;          % kern point: bot fiber side  [in]  (= r²/yb)
Fi   = F_des;            % initial prestress force [kip]
Fe   = eta * F_des;      % effective force after losses [kip]

rpt_file = fullfile(out_dir, 'Feasibility_Report.txt');
fid_r    = fopen(rpt_file, 'w');

w2('=================================================================\n', fid_r);
w2(' ECCENTRICITY BOUNDS — Kern-Point Linear Equations\n', fid_r);
w2(' Naaman Form:  e_o = kern +/- (1/F) * (M +/- sigma_allow * Z)\n', fid_r);
w2(sprintf(' F_i = %.2f kip  |  eta = %.2f  |  eta*F = %.2f kip\n', Fi, eta, Fe), fid_r);
w2('=================================================================\n\n', fid_r);

w2(sprintf('  Kern points:\n'), fid_r);
w2(sprintf('    k_t = St / Ac = %.2f / %.3f = %.4f in  (top-fiber side, r^2/yt)\n', St, Ac, kt_k), fid_r);
w2(sprintf('    k_b = Sb / Ac = %.2f / %.3f = %.4f in  (bot-fiber side, r^2/yb)\n', Sb, Ac, kb_k), fid_r);
w2(sprintf('    Z_t = St = %.2f in^3\n', St), fid_r);
w2(sprintf('    Z_b = Sb = %.2f in^3\n\n', Sb), fid_r);

w2(sprintf('  Allowable stresses:\n'), fid_r);
w2(sprintf('    Transfer: fci = +%.4f ksi (compr)  |  fti = +%.4f ksi (tens magnitude)\n', fci_allow, fti_allow), fid_r);
w2(sprintf('    Service:  fcs = +%.4f ksi (compr)  |  fts = +%.4f ksi (tens magnitude)\n\n', fcs_allow, fts_allow), fid_r);

w2('  General equations (Naaman Line labels in brackets):\n', fid_r);
w2('    TRANSFER (F_i):  -- Naaman notation, Table 4.2 --\n', fid_r);
w2('    [addl ]  e >= +k_t + (M_i - fci*Zt) / F_i      top compr <= +fci  [lower, additional]\n', fid_r);
w2('    [Line I] e <= +k_t + (M_i + fti*Zt) / F_i      top tens  >= -fti  [UPPER bound on e]\n', fid_r);
w2('    [addl ]  e >= -k_b + (M_i - fti*Zb) / F_i      bot tens  >= -fti  [lower, additional]\n', fid_r);
w2('    [LineII] e <= -k_b + (M_i + fci*Zb) / F_i      bot compr <= +fci  [UPPER bound on e]\n', fid_r);
w2('    SERVICE (eta*F):  -- Naaman notation, Table 4.2 --\n', fid_r);
w2('    [LineIV] e >= +k_t + (M_T - fcs*Zt) / (eta*F)  top compr <= +fcs  [LOWER bound on e]\n', fid_r);
w2('    [addl ]  e <= +k_t + (M_T + fts*Zt) / (eta*F)  top tens  >= -fts  [upper, additional]\n', fid_r);
w2('    [LineIII]e >= -k_b + (M_T - fts*Zb) / (eta*F)  bot tens  >= -fts  [LOWER bound on e]\n', fid_r);
w2('    [addl ]  e <= -k_b + (M_T + fcs*Zb) / (eta*F)  bot compr <= +fcs  [upper, additional]\n', fid_r);
w2('    Naaman governing:  e_min = max(Line IV, Line III)\n', fid_r);
w2('                       e_max = min(Line I,  Line II)\n\n', fid_r);

% Report at 4 key sections
x_key_pts = [0,  L/4,  L/2,  3*L/4];
x_lbl     = {'x = 0        (support)',              ...
             sprintf('x = L/4  = %.0f in', L/4),    ...
             sprintf('x = L/2  = %.0f in  (midspan)', L/2), ...
             sprintf('x = 3L/4 = %.0f in', 3*L/4)};

for ipt = 1:length(x_key_pts)
    xi   = x_key_pts(ipt);
    Mi_x = w_sw   * xi * (L - xi) / 2;
    MT_x = w_serv * xi * (L - xi) / 2;

    w2('------------------------------------------------------------------\n', fid_r);
    w2(sprintf('  AT  %s  (x/L = %.2f)\n', x_lbl{ipt}, xi/L), fid_r);
    w2('------------------------------------------------------------------\n', fid_r);
    w2(sprintf('  M_i = %.2f kip-in  |  M_T = %.2f kip-in\n\n', Mi_x, MT_x), fid_r);

    % Compute all 8 bounds
    eI    =  kt_k + (Mi_x - fci_allow*St) / Fi;
    eII   =  kt_k + (Mi_x + fti_allow*St) / Fi;
    eIII  = -kb_k + (Mi_x - fti_allow*Sb) / Fi;
    eIV   = -kb_k + (Mi_x + fci_allow*Sb) / Fi;
    eV    =  kt_k + (MT_x - fcs_allow*St) / Fe;
    eVI   =  kt_k + (MT_x + fts_allow*St) / Fe;
    eVII  = -kb_k + (MT_x - fts_allow*Sb) / Fe;
    eVIII = -kb_k + (MT_x + fcs_allow*Sb) / Fe;

    % ---- TRANSFER ---------------------------------------------------
    w2(sprintf('  TRANSFER  (F_i = %.2f kip):\n\n', Fi), fid_r);

    w2('  [addl]   Top fiber @ Transfer, compression (<=+fci):  e >= +k_t + (M_i - fci*Zt) / F_i   [lower bound, additional]\n', fid_r);
    w2(sprintf('         = +%.4f + (%.2f - %.4f x %.2f) / %.2f\n', kt_k, Mi_x, fci_allow, St, Fi), fid_r);
    w2(sprintf('         = +%.4f + (%.2f - %.2f) / %.2f\n', kt_k, Mi_x, fci_allow*St, Fi), fid_r);
    w2(sprintf('         = +%.4f + %.4f\n', kt_k, (Mi_x - fci_allow*St)/Fi), fid_r);
    w2(sprintf('         =  %.4f in\n\n', eI), fid_r);

    w2('  [Naaman Line I]   Top fiber @ Transfer, tension (>=-fti):  e <= +k_t + (M_i + fti*Zt) / F_i   [UPPER bound on e]\n', fid_r);
    w2(sprintf('         = +%.4f + (%.2f + %.4f x %.2f) / %.2f\n', kt_k, Mi_x, fti_allow, St, Fi), fid_r);
    w2(sprintf('         = +%.4f + (%.2f + %.2f) / %.2f\n', kt_k, Mi_x, fti_allow*St, Fi), fid_r);
    w2(sprintf('         = +%.4f + %.4f\n', kt_k, (Mi_x + fti_allow*St)/Fi), fid_r);
    w2(sprintf('         =  %.4f in\n\n', eII), fid_r);

    w2('  [addl]   Bot fiber @ Transfer, tension (>=-fti):  e >= -k_b + (M_i - fti*Zb) / F_i   [lower bound, additional]\n', fid_r);
    w2(sprintf('         = -%.4f + (%.2f - %.4f x %.2f) / %.2f\n', kb_k, Mi_x, fti_allow, Sb, Fi), fid_r);
    w2(sprintf('         = -%.4f + (%.2f - %.2f) / %.2f\n', kb_k, Mi_x, fti_allow*Sb, Fi), fid_r);
    w2(sprintf('         = -%.4f + %.4f\n', kb_k, (Mi_x - fti_allow*Sb)/Fi), fid_r);
    w2(sprintf('         =  %.4f in\n\n', eIII), fid_r);

    w2('  [Naaman Line II]  Bot fiber @ Transfer, compression (<=+fci):  e <= -k_b + (M_i + fci*Zb) / F_i   [UPPER bound on e]\n', fid_r);
    w2(sprintf('         = -%.4f + (%.2f + %.4f x %.2f) / %.2f\n', kb_k, Mi_x, fci_allow, Sb, Fi), fid_r);
    w2(sprintf('         = -%.4f + (%.2f + %.2f) / %.2f\n', kb_k, Mi_x, fci_allow*Sb, Fi), fid_r);
    w2(sprintf('         = -%.4f + %.4f\n', kb_k, (Mi_x + fci_allow*Sb)/Fi), fid_r);
    w2(sprintf('         =  %.4f in\n\n', eIV), fid_r);

    % ---- SERVICE ----------------------------------------------------
    w2(sprintf('  SERVICE  (eta*F = %.2f kip):\n\n', Fe), fid_r);

    w2('  [Naaman Line IV]  Top fiber @ Service, compression (<=+fcs):  e >= +k_t + (M_T - fcs*Zt) / (eta*F)  [LOWER bound on e]\n', fid_r);
    w2(sprintf('         = +%.4f + (%.2f - %.4f x %.2f) / %.2f\n', kt_k, MT_x, fcs_allow, St, Fe), fid_r);
    w2(sprintf('         = +%.4f + (%.2f - %.2f) / %.2f\n', kt_k, MT_x, fcs_allow*St, Fe), fid_r);
    w2(sprintf('         = +%.4f + %.4f\n', kt_k, (MT_x - fcs_allow*St)/Fe), fid_r);
    w2(sprintf('         =  %.4f in\n\n', eV), fid_r);

    w2('  [addl]   Top fiber @ Service, tension (>=-fts):  e <= +k_t + (M_T + fts*Zt) / (eta*F)  [upper bound, additional]\n', fid_r);
    w2(sprintf('         = +%.4f + (%.2f + %.4f x %.2f) / %.2f\n', kt_k, MT_x, fts_allow, St, Fe), fid_r);
    w2(sprintf('         = +%.4f + (%.2f + %.2f) / %.2f\n', kt_k, MT_x, fts_allow*St, Fe), fid_r);
    w2(sprintf('         = +%.4f + %.4f\n', kt_k, (MT_x + fts_allow*St)/Fe), fid_r);
    w2(sprintf('         =  %.4f in\n\n', eVI), fid_r);

    w2('  [Naaman Line III] Bot fiber @ Service, tension (>=-fts):  e >= -k_b + (M_T - fts*Zb) / (eta*F)  [LOWER bound on e]\n', fid_r);
    w2(sprintf('         = -%.4f + (%.2f - %.4f x %.2f) / %.2f\n', kb_k, MT_x, fts_allow, Sb, Fe), fid_r);
    w2(sprintf('         = -%.4f + (%.2f - %.2f) / %.2f\n', kb_k, MT_x, fts_allow*Sb, Fe), fid_r);
    w2(sprintf('         = -%.4f + %.4f\n', kb_k, (MT_x - fts_allow*Sb)/Fe), fid_r);
    w2(sprintf('         =  %.4f in\n\n', eVII), fid_r);

    w2('  [addl]   Bot fiber @ Service, compression (<=+fcs):  e <= -k_b + (M_T + fcs*Zb) / (eta*F)  [upper bound, additional]\n', fid_r);
    w2(sprintf('         = -%.4f + (%.2f + %.4f x %.2f) / %.2f\n', kb_k, MT_x, fcs_allow, Sb, Fe), fid_r);
    w2(sprintf('         = -%.4f + (%.2f + %.2f) / %.2f\n', kb_k, MT_x, fcs_allow*Sb, Fe), fid_r);
    w2(sprintf('         = -%.4f + %.4f\n', kb_k, (MT_x + fcs_allow*Sb)/Fe), fid_r);
    w2(sprintf('         =  %.4f in\n\n', eVIII), fid_r);

    % ---- Governing bounds summary -----------------------------------
    % Lower bounds on e (e must be >= these values):
    %   addl  = eI  (transfer, top compr)
    %   LineIV= eV  (service, top compr)   [Naaman Line IV]
    %   addl  = eIII(transfer, bot tens)
    %   LineIII=eVII(service, bot tens)    [Naaman Line III]
    % Upper bounds on e (e must be <= these values):
    %   LineI = eII (transfer, top tens)   [Naaman Line I]
    %   LineII= eIV (transfer, bot compr)  [Naaman Line II]
    %   addl  = eVI (service, top tens)
    %   addl  = eVIII(service, bot compr)
    e_lowers = [eI,   eV,   eIII, eVII];
    e_uppers = [eII,  eIV,  eVI,  eVIII];
    l_names  = {'addl (xfer top compr)', 'Naaman Line IV (svc top compr)', ...
                'addl (xfer bot tens)',  'Naaman Line III (svc bot tens)'};
    u_names  = {'Naaman Line I (xfer top tens)',  'Naaman Line II (xfer bot compr)', ...
                'addl (svc top tens)',            'addl (svc bot compr)'};
    [e_min_gov, i_lo] = max(e_lowers);
    [e_max_gov, i_hi] = min(e_uppers);

    w2('  GOVERNING BOUNDS (Naaman Lines I-IV + additional):\n', fid_r);
    w2(sprintf('    Lower bounds on e:\n'), fid_r);
    w2(sprintf('      [addl, xfer top compr] eI   = %.4f in\n', eI), fid_r);
    w2(sprintf('      [Naaman Line IV]       eV   = %.4f in  (service, top compr)\n', eV), fid_r);
    w2(sprintf('      [addl, xfer bot tens]  eIII = %.4f in\n', eIII), fid_r);
    w2(sprintf('      [Naaman Line III]      eVII = %.4f in  (service, bot tens)\n', eVII), fid_r);
    w2(sprintf('    --> e_min = max of above = %.4f in   [%s governs]\n', e_min_gov, l_names{i_lo}), fid_r);
    w2(sprintf('    Upper bounds on e:\n'), fid_r);
    w2(sprintf('      [Naaman Line I]        eII  = %.4f in  (transfer, top tens)\n', eII), fid_r);
    w2(sprintf('      [Naaman Line II]       eIV  = %.4f in  (transfer, bot compr)\n', eIV), fid_r);
    w2(sprintf('      [addl, svc top tens]   eVI  = %.4f in\n', eVI), fid_r);
    w2(sprintf('      [addl, svc bot compr]  eVIII= %.4f in\n', eVIII), fid_r);
    w2(sprintf('    --> e_max = min of above = %.4f in   [%s governs]\n', e_max_gov, u_names{i_hi}), fid_r);

    if e_max_gov > e_min_gov
        w2(sprintf('    Feasible range: %.4f <= e <= %.4f in  (width = %.4f in)\n', ...
                   e_min_gov, e_max_gov, e_max_gov - e_min_gov), fid_r);
    else
        w2(sprintf('    *** NO FEASIBLE RANGE at this section (F = %.0f kip) ***\n', Fi), fid_r);
    end

    if abs(xi - L/2) < 1   % midspan: check design point
        w2(sprintf('\n    Design point check (midspan): e_des = %.4f in\n', e_des), fid_r);
        if e_des >= e_min_gov - 1e-6 && e_des <= e_max_gov + 1e-6
            w2('    --> e_des is INSIDE the feasible range   OK\n\n', fid_r);
        else
            w2('    --> e_des is OUTSIDE the feasible range  FAIL\n\n', fid_r);
        end
    else
        w2('\n', fid_r);
    end
end

fclose(fid_r);
fprintf('\n  Eccentricity bounds report saved to:\n    %s\n\n', rpt_file);

%% =========================================================================
%% 8.  SAVE FIGURES  →  Project1/
%% =========================================================================
figs  = [fig1, fig2];
names = {'Figure_1_Magnel_Diagram', 'Figure_2_Feasibility_Zone'};
for k = 1:numel(figs)
    base = fullfile(out_dir, names{k});
    saveas(figs(k), [base '.png']);
    savefig(figs(k), [base '.fig']);
    fprintf('  Saved: %s (.png + .fig)\n', names{k});
end

fprintf('\n=================================================================\n');
fprintf(' Feasibility Design Chart complete.  Two figures saved to:\n');
fprintf('   %s\n', out_dir);
fprintf('=================================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function [A, I, yc, yt, yb] = computeSectionProps(v)
% Compute area, moment of inertia, and centroid from polygon vertices
% using the shoelace (Green's theorem) formula.
% Vertices must be ordered counter-clockwise.
n = size(v, 1);
v = [v; v(1,:)];      % close polygon

% Area
A = 0;
for i = 1:n
    A = A + (v(i,1)*v(i+1,2) - v(i+1,1)*v(i,2));
end
A = abs(A) / 2;

% Centroid
Cx = 0;  Cy = 0;
for i = 1:n
    f  = v(i,1)*v(i+1,2) - v(i+1,1)*v(i,2);
    Cx = Cx + (v(i,1) + v(i+1,1)) * f;
    Cy = Cy + (v(i,2) + v(i+1,2)) * f;
end
Cx = Cx / (6*A);
Cy = Cy / (6*A);

% Moment of inertia about centroidal x-axis
Ix = 0;
for i = 1:n
    xi  = v(i,1)   - Cx;   yi  = v(i,2)   - Cy;
    xi1 = v(i+1,1) - Cx;   yi1 = v(i+1,2) - Cy;
    f   = xi*yi1 - xi1*yi;
    Ix  = Ix + (yi^2 + yi*yi1 + yi1^2) * f;
end
I  = abs(Ix) / 12;

yc = Cy;
yt = max(v(:,2)) - Cy;    % centroid to top  fiber
yb = Cy - min(v(:,2));    % centroid to bottom fiber
end

% -------------------------------------------------------------------------
function w2(str, fid)
% Write str to both the console and an open file fid.
% Use fprintf(str) so that \n in literal strings is treated as newline.
fprintf(str);
fprintf(fid, str);
end

% -------------------------------------------------------------------------
function s = chkS(f, fmin, fmax)
if f >= fmin - 1e-8 && f <= fmax + 1e-8
    s = 'OK';
else
    s = '*** FAIL ***';
end
end
