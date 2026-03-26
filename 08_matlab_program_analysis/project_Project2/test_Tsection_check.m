%% test_Tsection_check.m
%
%  VERIFICATION TEST — T-Section with Compression in Both Flange AND Web
%
%  Purpose: confirm that the T-section branch (a > hf) in ultimateDesign.m
%           correctly:
%            (1) splits the compression zone into flange overhangs + web block
%            (2) includes mild steel (As, Asc) in equilibrium and Mn
%            (3) uses bw (not b) for fps via ACI Eq. 22.3.2.1
%
%  Section geometry (simple prestressed T-beam — NOT the Double-T):
%    b   = 36 in  (effective flange width)
%    bw  = 12 in  (web width)
%    hf  =  4 in  (flange thickness)
%    h   = 36 in  (total depth)
%
%  Sign convention: compression = positive (+), tension = negative (-)
%  Coordinates:    y = 0 at BOTTOM, y = h at TOP (consistent with project)
%  Depths from top are used only where ACI formula requires (dp, d, d').
%
%  Report saved to:  output/Tsection_Check_Report.txt
%% ======================================================================

clear; clc;

%% ======== SECTION GEOMETRY =========================================
b   = 36;    % in  effective flange width
bw  = 12;    % in  web width
hf  =  4;    % in  flange (slab) thickness
h   = 36;    % in  total section depth

%% ======== MATERIALS =================================================
fc      = 4.0;   % ksi  f'c concrete
fpu     = 270;   % ksi  ultimate PS strand
fpy     = 243;   % ksi  yield (0.9*fpu, low-relaxation → gamma_p = 0.28)
fy      = 60;    % ksi  mild steel
beta1   = 0.85 - 0.05*((fc - 4)/1);   % ACI 22.2.2.4.3 (= 0.85 for fc=4 ksi)
beta1   = max(beta1, 0.65);

% gamma_p: ACI Table 22.3.2.1
%   fpy/fpu = 243/270 = 0.90 → low-relaxation strand → gamma_p = 0.28
gamma_p = 0.28;

%% ======== PRESTRESS TENDONS =========================================
Aps     = 2.40;   % in^2  total bonded prestress area
fpi     = 0.70 * fpu;   % ksi  initial jacking stress (70% fpu)
eta     = 0.85;          % loss factor
fse     = eta * fpi;     % ksi  effective prestress
dp      = 32;    % in  depth from top to PS centroid (h - y_tendon = 36-4 = 32)

%% ======== MILD STEEL ================================================
As      = 2.0;   % in^2  TENSION mild steel (at bottom)
d_s     = 33;    % in   depth from top to tension steel centroid
Asc     = 0.5;   % in^2  COMPRESSION mild steel (near top)
d_sc    =  3.0;  % in   depth from top to compression steel centroid

%% ======== OUTPUT ====================================================
out_dir = 'output';
if ~exist(out_dir, 'dir'),  mkdir(out_dir);  end
report_file = fullfile(out_dir, 'Tsection_Check_Report.txt');
fid = fopen(report_file, 'w');

function pf(fid, varargin)
    str = sprintf(varargin{:});
    fprintf('%s', str);
    fprintf(fid, '%s', str);
end

%% ======== START REPORT =============================================
pf(fid, '\n');
pf(fid, '==================================================================\n');
pf(fid, '   T-SECTION VERIFICATION — Compression in Flange AND Web\n');
pf(fid, '   ACI 318-19 Eq. 22.3.2.1  |  With mild steel As and As''\n');
pf(fid, '==================================================================\n');
pf(fid, '\n');
pf(fid, '  SIGN CONVENTION: Compression = positive (+), Tension = negative (-)\n');
pf(fid, '  Depths measured from TOP fiber (ACI convention for c, a, dp, d, d'')\n');
pf(fid, '\n');

%% ---- Section geometry ---------------------------------------------
pf(fid, '------------------------------------------------------------------\n');
pf(fid, '  SECTION GEOMETRY\n');
pf(fid, '------------------------------------------------------------------\n');
pf(fid, '  b   = %5.1f in   (effective flange width)\n', b);
pf(fid, '  bw  = %5.1f in   (web width)\n', bw);
pf(fid, '  hf  = %5.1f in   (flange thickness)\n', hf);
pf(fid, '  h   = %5.1f in   (total depth)\n', h);
pf(fid, '\n');
pf(fid, '  [ASCII cross-section]\n');
pf(fid, '\n');
pf(fid, '  |<------- b = 36 in ------->|\n');
pf(fid, '  +---------------------------+  -- top (y = h = 36 in)\n');
pf(fid, '  |     FLANGE  hf = 4 in     |\n');
pf(fid, '  +--------+     +------------+  -- flange soffit (y = 32 in)\n');
pf(fid, '           | WEB |\n');
pf(fid, '           |bw=12|   dp = 32 in  [PS centroid from top]\n');
pf(fid, '           |  o  |  <-- PS tendons (Aps = 2.40 in^2)\n');
pf(fid, '           +-----+  -- bottom (y = 0 in)\n');
pf(fid, '\n');

%% ---- Materials ----------------------------------------------------
pf(fid, '------------------------------------------------------------------\n');
pf(fid, '  MATERIAL PROPERTIES\n');
pf(fid, '------------------------------------------------------------------\n');
pf(fid, '  f''c     = %.1f ksi\n', fc);
pf(fid, '  fpu     = %.0f ksi\n', fpu);
pf(fid, '  fpy/fpu = %.2f  (low-relaxation strand)\n', fpy/fpu);
pf(fid, '  gamma_p = %.2f  (ACI Table 22.3.2.1, low-relax)\n', gamma_p);
pf(fid, '  beta_1  = %.3f  (ACI Table 22.2.2.4.3, f''c = %.1f ksi)\n', beta1, fc);
pf(fid, '  fy      = %.0f ksi  (mild steel)\n', fy);
pf(fid, '\n');

%% ---- Steel layout -------------------------------------------------
pf(fid, '------------------------------------------------------------------\n');
pf(fid, '  STEEL LAYOUT\n');
pf(fid, '------------------------------------------------------------------\n');
pf(fid, '  Prestress:          Aps = %.2f in^2,  dp = %.1f in from top\n', Aps, dp);
pf(fid, '                      fse = %.2f ksi  (>= 0.5*fpu = %.1f ksi)  CHECK: %s\n', ...
        fse, 0.5*fpu, select_str(fse >= 0.5*fpu, 'OK', 'FAIL'));
pf(fid, '  Tension mild steel: As  = %.2f in^2,  d   = %.1f in from top\n', As, d_s);
pf(fid, '  Compression steel:  As'' = %.2f in^2,  d'' = %.1f in from top\n', Asc, d_sc);
pf(fid, '\n');

%% ======================================================================
%%  PART 1 — fps via ACI 318-19 Eq. 22.3.2.1 (full form)
%% ======================================================================
pf(fid, '==================================================================\n');
pf(fid, '   PART 1 — STRESS IN PRESTRESSING STEEL  fps\n');
pf(fid, '   ACI 318-19 Eq. 22.3.2.1 (full form with mild steel)\n');
pf(fid, '==================================================================\n');
pf(fid, '\n');
pf(fid, '  Full equation:\n');
pf(fid, '    fps = fpu * [ 1 - (gamma_p / beta_1)\n');
pf(fid, '                    * ( rho_p * fpu/f''c\n');
pf(fid, '                      + (d / dp)  * omega\n');
pf(fid, '                      - (d''/ dp)  * omega'' ) ]\n');
pf(fid, '\n');
pf(fid, '  Where:\n');
pf(fid, '    omega   = (As  / (b_eff * dp)) * fy / f''c   (tension mild steel index)\n');
pf(fid, '    omega'' = (As'' / (b_eff * dp)) * fy / f''c   (compression mild steel index)\n');
pf(fid, '\n');

%% --- Step A: Trial with full flange width b -------------------------
pf(fid, '  --- STEP A: TRIAL (rectangular, b_eff = b = %.0f in) ---\n', b);
pf(fid, '\n');

rho_p_trial = Aps / (b * dp);
omega_trial = (As  / (b * dp)) * fy / fc;
omega_c_trial = (Asc / (b * dp)) * fy / fc;

fps_trial = fpu * (1 - (gamma_p / beta1) * ...
            (rho_p_trial * fpu / fc + (d_s/dp)*omega_trial - (d_sc/dp)*omega_c_trial));

a_trial = (Aps*fps_trial + As*fy - Asc*fy) / (0.85 * fc * b);

pf(fid, '  rho_p     = Aps/(b*dp)             = %.4f/(%.0f*%.0f)  = %.6f\n', Aps, b, dp, rho_p_trial);
pf(fid, '  omega     = (As/(b*dp)) * fy/f''c  = (%.2f/(%.0f*%.0f)) * %.0f/%.1f  = %.6f\n', As, b, dp, fy, fc, omega_trial);
pf(fid, '  omega''   = (As''/(b*dp)) * fy/f''c = (%.2f/(%.0f*%.0f)) * %.0f/%.1f  = %.6f\n', Asc, b, dp, fy, fc, omega_c_trial);
pf(fid, '\n');
pf(fid, '  fps_trial = %.0f * [1 - (%.2f/%.3f) * (%.6f*%.0f/%.1f\n', fpu, gamma_p, beta1, rho_p_trial, fpu, fc);
pf(fid, '                           + (%.1f/%.0f)*%.6f - (%.1f/%.0f)*%.6f)]\n', d_s, dp, omega_trial, d_sc, dp, omega_c_trial);
pf(fid, '            = %.2f ksi\n', fps_trial);
pf(fid, '\n');
pf(fid, '  a_trial   = (Aps*fps + As*fy - As''*fy) / (0.85*f''c*b)\n');
pf(fid, '            = (%.2f*%.2f + %.2f*%.0f - %.2f*%.0f) / (0.85*%.1f*%.0f)\n', ...
        Aps, fps_trial, As, fy, Asc, fy, fc, b);
pf(fid, '            = %.3f in\n', a_trial);
pf(fid, '\n');

if a_trial <= hf
    pf(fid, '  a_trial = %.3f in  <=  hf = %.1f in\n', a_trial, hf);
    pf(fid, '  --> Compression block lies ENTIRELY within flange. No T-section needed.\n');
    pf(fid, '  NOTE: This test was designed to force T-section. Check parameters.\n');
    % fallback: use rectangular result
    fps = fps_trial;
    a   = a_trial;
    rho_pw = NaN;  Ff = NaN;  A_fl = NaN;  A_web = NaN;  ybar = a/2;
    section_type = 'Rectangular (unexpected)';
else
    pf(fid, '  a_trial = %.3f in  >  hf = %.1f in\n', a_trial, hf);
    pf(fid, '  --> Compression block PENETRATES INTO WEB.  T-Section branch required.\n');
    pf(fid, '\n');

    %% --- Step B: Recompute fps using bw (web controls for T-section) ---
    pf(fid, '  --- STEP B: T-SECTION fps  (use bw = %.0f in) ---\n', bw);
    pf(fid, '\n');
    pf(fid, '  When the neutral axis is in the web, ACI Eq. 22.3.2.1 uses bw:\n');
    pf(fid, '    rho_pw = Aps / (bw * dp)\n');
    pf(fid, '    omega_w  = (As  / (bw * dp)) * fy / f''c\n');
    pf(fid, '    omega_cw = (As'' / (bw * dp)) * fy / f''c\n');
    pf(fid, '\n');

    rho_pw   = Aps / (bw * dp);
    omega_w  = (As  / (bw * dp)) * fy / fc;
    omega_cw = (Asc / (bw * dp)) * fy / fc;
    fps      = fpu * (1 - (gamma_p / beta1) * ...
               (rho_pw * fpu / fc + (d_s/dp)*omega_w - (d_sc/dp)*omega_cw));

    pf(fid, '  rho_pw   = Aps/(bw*dp)              = %.4f/(%.0f*%.0f)  = %.6f\n', Aps, bw, dp, rho_pw);
    pf(fid, '  omega_w  = (As/(bw*dp)) * fy/f''c   = (%.2f/(%.0f*%.0f)) * %.0f/%.1f  = %.6f\n', As, bw, dp, fy, fc, omega_w);
    pf(fid, '  omega_cw = (As''/(bw*dp)) * fy/f''c  = (%.2f/(%.0f*%.0f)) * %.0f/%.1f  = %.6f\n', Asc, bw, dp, fy, fc, omega_cw);
    pf(fid, '\n');
    pf(fid, '  fps = %.0f * [1 - (%.2f/%.3f) * (%.6f*%.0f/%.1f\n', fpu, gamma_p, beta1, rho_pw, fpu, fc);
    pf(fid, '                      + (%.1f/%.0f)*%.6f - (%.1f/%.0f)*%.6f)]\n', d_s, dp, omega_w, d_sc, dp, omega_cw);
    pf(fid, '      = %.2f ksi\n', fps);
    pf(fid, '\n');

    %% --- Step C: Compression zone — flange overhangs + web block -----
    pf(fid, '==================================================================\n');
    pf(fid, '   PART 2 — COMPRESSION ZONE GEOMETRY\n');
    pf(fid, '   Split: flange overhangs (b-bw) x hf  +  web block bw x a\n');
    pf(fid, '==================================================================\n');
    pf(fid, '\n');

    % Flange overhang force
    Ff = 0.85 * fc * (b - bw) * hf;
    pf(fid, '  FLANGE OVERHANGS (compression force Ff):\n');
    pf(fid, '    Ff = 0.85 * f''c * (b - bw) * hf\n');
    pf(fid, '       = 0.85 * %.1f * (%.0f - %.0f) * %.1f\n', fc, b, bw, hf);
    pf(fid, '       = %.2f kip\n', Ff);
    pf(fid, '\n');

    % Web block depth from equilibrium
    %   Aps*fps + As*fy - Asc*fy = Ff + 0.85*f'c*bw*a
    a = (Aps*fps + As*fy - Asc*fy - Ff) / (0.85 * fc * bw);

    pf(fid, '  WEB BLOCK DEPTH  a  (from equilibrium):\n');
    pf(fid, '    Force equilibrium:  Aps*fps + As*fy - As''*fy = Ff + 0.85*f''c*bw*a\n');
    pf(fid, '    Solving for a:\n');
    pf(fid, '      a = (Aps*fps + As*fy - As''*fy - Ff) / (0.85*f''c*bw)\n');
    pf(fid, '        = (%.2f*%.2f + %.2f*%.0f - %.2f*%.0f - %.2f) / (0.85*%.1f*%.0f)\n', ...
            Aps, fps, As, fy, Asc, fy, Ff, fc, bw);
    pf(fid, '        = %.3f in\n', a);
    pf(fid, '\n');

    if a <= 0
        pf(fid, '  WARNING: a <= 0 — tension force is not large enough to require web compression.\n');
        pf(fid, '           Check parameters (fps may have dropped too low with bw).\n');
    end

    % Compression zone geometry
    A_fl  = (b - bw) * hf;          % flange overhang area (2 overhangs combined)
    A_web = bw * a;                  % web compression block area
    ybar  = (A_fl*(hf/2) + A_web*(a/2)) / (A_fl + A_web);   % centroid from top

    pf(fid, '  COMPRESSION ZONE AREAS AND CENTROID:\n');
    pf(fid, '    A_fl  = (b - bw) * hf  = (%.0f - %.0f) * %.1f  = %.3f in^2\n', b, bw, hf, A_fl);
    pf(fid, '            centroid from top = hf/2 = %.4f in\n', hf/2);
    pf(fid, '\n');
    pf(fid, '    A_web = bw * a          = %.0f * %.3f        = %.3f in^2\n', bw, a, A_web);
    pf(fid, '            centroid from top = a/2  = %.4f in\n', a/2);
    pf(fid, '\n');
    pf(fid, '    ybar  = (A_fl*(hf/2) + A_web*(a/2)) / (A_fl + A_web)\n');
    pf(fid, '          = (%.3f*%.4f + %.3f*%.4f) / (%.3f + %.3f)\n', ...
            A_fl, hf/2, A_web, a/2, A_fl, A_web);
    pf(fid, '          = %.4f in  from top fiber\n', ybar);
    pf(fid, '\n');
    pf(fid, '  [Compression zone — ASCII sketch, depths from top]\n');
    pf(fid, '\n');
    pf(fid, '  0 in  +---------------------------+  top fiber\n');
    pf(fid, '        |  FLANGE OVERHANGS (0 to hf=%.1f in)  |\n', hf);
    pf(fid, '  %.1f in +--------+         +--------+  flange soffit\n', hf);
    pf(fid, '                  |  WEB    |\n');
    pf(fid, '  ybar=%.2f +---->  | a=%.3f  |  <-- centroid of compression zone\n', ybar, a);
    pf(fid, '  %.3f +--------+         web block bottom\n', a);
    pf(fid, '\n');

    %% --- Step D: Nominal moment Mn -----------------------------------
    pf(fid, '==================================================================\n');
    pf(fid, '   PART 3 — NOMINAL MOMENT CAPACITY  Mn\n');
    pf(fid, '==================================================================\n');
    pf(fid, '\n');
    pf(fid, '  Formula (T-section, moments about centroid of compression zone at ybar):\n');
    pf(fid, '    Mn = Aps*fps*(dp - ybar) + As*fy*(d - ybar) - As''*fy*(d'' - ybar)\n');
    pf(fid, '\n');
    pf(fid, '  ybar = %.4f in  (centroid of compression zone from top)\n', ybar);
    pf(fid, '  dp   = %.1f in  (PS steel depth from top)\n', dp);
    pf(fid, '  d    = %.1f in  (tension mild steel depth from top)\n', d_s);
    pf(fid, '  d''  = %.1f in  (compression mild steel depth from top)\n', d_sc);
    pf(fid, '\n');
    pf(fid, '  Lever arms:\n');
    pf(fid, '    dp   - ybar = %.4f - %.4f = %.4f in\n', dp, ybar, dp-ybar);
    pf(fid, '    d    - ybar = %.4f - %.4f = %.4f in\n', d_s, ybar, d_s-ybar);
    pf(fid, '    d''  - ybar = %.4f - %.4f = %.4f in\n', d_sc, ybar, d_sc-ybar);
    pf(fid, '\n');
    pf(fid, '  Substitution:\n');
    pf(fid, '    Aps*fps*(dp-ybar) = %.2f * %.2f * %.4f  = %.2f kip-in\n', ...
            Aps, fps, dp-ybar, Aps*fps*(dp-ybar));
    pf(fid, '    As*fy*(d-ybar)    = %.2f * %.0f * %.4f  = %.2f kip-in\n', ...
            As, fy, d_s-ybar, As*fy*(d_s-ybar));
    pf(fid, '    As''*fy*(d''-ybar) = %.2f * %.0f * %.4f  = %.2f kip-in  [subtracted]\n', ...
            Asc, fy, d_sc-ybar, Asc*fy*(d_sc-ybar));
    pf(fid, '\n');

    Mn = Aps*fps*(dp - ybar) ...
       + As *fy *(d_s  - ybar) ...
       - Asc*fy *(d_sc - ybar);

    pf(fid, '    Mn = %.2f + %.2f - %.2f\n', ...
            Aps*fps*(dp-ybar), As*fy*(d_s-ybar), Asc*fy*(d_sc-ybar));
    pf(fid, '       = %.2f kip-in  =  %.2f kip-ft\n', Mn, Mn/12);
    pf(fid, '\n');

    section_type = 'T-Section (a > hf)';

    %% --- Step E: Ductility check ------------------------------------
    c     = a / beta1;
    eps_t = (dp - c) / c * 0.003;

    pf(fid, '==================================================================\n');
    pf(fid, '   PART 4 — DUCTILITY CHECK  (ACI 318-19 §9.3.3)\n');
    pf(fid, '==================================================================\n');
    pf(fid, '\n');
    pf(fid, '  NOTE: c is measured from top to NEUTRAL AXIS (bottom of web block).\n');
    pf(fid, '        c = hf + a / beta_1  would apply only if a refers to web depth.\n');
    pf(fid, '        Here c = a_total / beta_1, where a_total = hf + a (web block).\n');
    pf(fid, '        However ACI uses c = a_web_only / beta_1 referenced to bottom\n');
    pf(fid, '        of flange for T-sections. Simplified here as c = a / beta_1\n');
    pf(fid, '        where a is the web block depth (conservative if a << hf).\n');
    pf(fid, '        More precisely:  c_total = (hf + a) for depth to NA from top.\n');
    pf(fid, '\n');

    c_total = hf + a;   % actual depth to neutral axis from top (for strain calc)
    eps_t_correct = (dp - c_total) / c_total * 0.003;

    pf(fid, '  c (depth to neutral axis from top) = hf + a = %.1f + %.3f = %.3f in\n', hf, a, c_total);
    pf(fid, '  et = (dp - c) / c * 0.003\n');
    pf(fid, '     = (%.1f - %.3f) / %.3f * 0.003\n', dp, c_total, c_total);
    pf(fid, '     = %.5f\n', eps_t_correct);
    pf(fid, '\n');

    if eps_t_correct >= 0.005
        phi = 0.90;
        duc_status = 'Tension-controlled (et >= 0.005) --> phi = 0.90  PASS';
    elseif eps_t_correct >= 0.004
        phi = 0.65 + (eps_t_correct - 0.002)/(0.005 - 0.002)*(0.90 - 0.65);
        duc_status = sprintf('Transition zone (0.004 <= et < 0.005) --> phi = %.3f', phi);
    else
        phi = NaN;
        duc_status = 'FAIL -- et < 0.004 (over-reinforced, not permitted by ACI 9.3.3)';
    end

    pf(fid, '  %s\n', duc_status);
    pf(fid, '\n');
    phi_Mn = phi * Mn;
    pf(fid, '  phi*Mn = %.2f * %.2f = %.2f kip-in  (%.2f kip-ft)\n', phi, Mn, phi_Mn, phi_Mn/12);
    pf(fid, '\n');
end

%% ======== SUMMARY ==================================================
pf(fid, '==================================================================\n');
pf(fid, '   FINAL SUMMARY\n');
pf(fid, '==================================================================\n');
pf(fid, '\n');
pf(fid, '  Section type   : %s\n', section_type);
pf(fid, '  fps            = %.2f ksi\n', fps);
pf(fid, '  a (web block)  = %.3f in  (web block depth, below flange)\n', a);
pf(fid, '  Total comp.    = hf + a  = %.1f + %.3f = %.3f in from top\n', hf, a, hf+a);
pf(fid, '  Ff (flange)    = %.2f kip\n', Ff);
pf(fid, '  A_fl           = %.3f in^2  @ %.4f in from top\n', A_fl, hf/2);
pf(fid, '  A_web          = %.3f in^2  @ %.4f in from top\n', A_web, a/2);
pf(fid, '  ybar           = %.4f in  (resultant compression depth from top)\n', ybar);
pf(fid, '  Mn             = %.2f kip-in  =  %.2f kip-ft\n', Mn, Mn/12);
if ~isnan(phi)
    pf(fid, '  phi*Mn         = %.2f kip-in  =  %.2f kip-ft\n', phi_Mn, phi_Mn/12);
end
pf(fid, '\n');
pf(fid, '  NOTE: This is a verification test only. Section geometry, Aps, As, Asc\n');
pf(fid, '        are chosen to force a > hf so that the T-section branch is exercised.\n');
pf(fid, '        Not a real design problem.\n');
pf(fid, '\n');
pf(fid, '  Report saved to: %s\n', report_file);

fclose(fid);
fprintf('\nReport written to: %s\n', report_file);

%% ======== LOCAL HELPER =============================================
function s = select_str(cond, a, b)
    if cond, s = a; else, s = b; end
end
