%% RUNALLCASES.m
%  Runs feasibility4Sections.m against three beam configurations to
%  demonstrate all three feasibility outcomes:
%
%    Case 1 — Normal          : optimized F found, kern-point analysis runs
%    Case 2 — High Live Load  : no zone at e_cg (service controls dominate)
%    Case 3 — Tendon Too High : no zone at e_cg (bad tendon placement)
%
%  Output figures and reports are saved to separate subfolders:
%    output/Case1_Normal/
%    output/Case2_HighLoad_Infeasible/
%    output/Case3_BadTendon_NoZone/
%
%  USAGE: Run this script from MATLAB with working directory set to
%         the folder containing feasibility4Sections.m (project_Project3 Feasibility/).

clear;

% script_dir is shared with feasibility4Sections.m (both in same folder).
% Setting it here prevents mfilename from returning the wrong path
% when feasibility4Sections.m is called via run().
script_dir = fileparts(mfilename('fullpath'));

fs_script = fullfile(script_dir, 'feasibility4Sections.m');

%% =========================================================================
%% CASE 1 — Normal Design  (Project 3 as-is)
%%   Expected: F_req = 424.5 kip  (governing from midspan, Line IV)
%%             Global feasibility: OK (marginal at drape point)
%%             Kern-point analysis: runs for all 4 sections
%% =========================================================================
fprintf('\n');
fprintf('###############################################################\n');
fprintf('##  CASE 1 — Normal Design (Project 3, optimized F = 424.5 kip)\n');
fprintf('###############################################################\n\n');

case_tag      = 'Case1_Normal';
w_LL_override = [];      % use inputData value (0.035 kip/in)
tendon_y_all  = [];      % use inputData tendon profile

run(fs_script);
close all;

%% =========================================================================
%% CASE 2 — Very High Live Load  (5x normal LL)
%%   w_LL = 5 x 0.035 = 0.175 kip/in  (= 2100 lb/ft, vs normal 420 lb/ft)
%%
%%   Expected: lineIII (service bottom tension upper bound) drops BELOW
%%             lineII  (transfer compression lower bound) at most sections.
%%             => Fmin = Inf at critical sections
%%             => "NO OPTIMIZED FORCE FOUND" — kern-point skipped.
%%
%%   Physical interpretation: service moment demand so large that no
%%   prestress force can simultaneously satisfy transfer compression AND
%%   service tension with the given tendon layout.
%% =========================================================================
fprintf('\n');
fprintf('###############################################################\n');
fprintf('##  CASE 2 — 5x Live Load (No Feasible F — Service Dominates)\n');
fprintf('###############################################################\n\n');

case_tag      = 'Case2_HighLoad_Infeasible';
w_LL_override = 5 * (0.42 / 12);    % 5 x 0.035 = 0.175 kip/in
tendon_y_all  = [];

run(fs_script);
close all;

%% =========================================================================
%% CASE 3 — Tendon Placed Above Centroid  (all tendons at y = 26 in)
%%   y_tendon = 26 in  =>  e_cg = yc - 26 = 20.362 - 26 = -5.638 in
%%   (tendon group is ABOVE the centroid — physically incorrect for a
%%    simply-supported prestressed beam)
%%
%%   Expected: at midspan and intermediate sections, the Magnel feasible
%%             zone exists at positive eccentricities but the tendon is at
%%             e_cg = -5.638 in (far left of the feasible zone).
%%             => Fmin = Inf at sections 2, 3, 4
%%             => "NO OPTIMIZED FORCE FOUND" — kern-point skipped.
%%
%%   Physical interpretation: tendon above centroid induces upward
%%   prestress eccentricity, adding tension to the bottom fiber at
%%   transfer rather than compression — incompatible with service demand.
%% =========================================================================
fprintf('\n');
fprintf('###############################################################\n');
fprintf('##  CASE 3 — Tendon Above Centroid  (e_cg = -5.64 in, No Zone)\n');
fprintf('###############################################################\n\n');

case_tag      = 'Case3_BadTendon_NoZone';
w_LL_override = [];
tendon_y_all  = 26;     % all 4 tendons at y=26 in (flange soffit, above centroid)

run(fs_script);
close all;

%% =========================================================================
%% DONE
%% =========================================================================
fprintf('\n');
fprintf('###############################################################\n');
fprintf('##  ALL CASES COMPLETE\n');
fprintf('##\n');
fprintf('##  Results saved to:\n');
fprintf('##    output/Case1_Normal/\n');
fprintf('##    output/Case2_HighLoad_Infeasible/\n');
fprintf('##    output/Case3_BadTendon_NoZone/\n');
fprintf('###############################################################\n\n');
