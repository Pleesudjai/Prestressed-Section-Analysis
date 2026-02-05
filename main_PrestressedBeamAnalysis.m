%% MAIN_PRESTRESSED_BEAM_ANALYSIS
% Main driver script for prestressed concrete beam analysis
% Creates N, V, M diagrams with detailed visualization
%
% Author: Generated for structural analysis
% Date: 2025
%
% Usage:
%   1. Edit inputPrestressedBeam.m to define your beam
%   2. Run this script
%   3. Results and plots will be generated automatically

clear; clc; close all;

%% Load input data
fprintf('Loading input data...\n');
[beam, section, materials, prestress, reinforcement, loads] = inputPrestressedBeam_HW1();

%% Run analysis
fprintf('\nRunning analysis...\n');
results = analyzePrestressedBeam(beam, section, materials, prestress, reinforcement, loads);

%% Display key results
fprintf('\n========================================\n');
fprintf('        ANALYSIS RESULTS SUMMARY\n');
fprintf('========================================\n');

% Maximum forces
[M_max, idx_M_max] = max(abs(results.M));
[V_max, idx_V_max] = max(abs(results.V));
[N_max, ~] = max(abs(results.N));

fprintf('\nMaximum Internal Forces:\n');
fprintf('  M_max = %.0f kip-in (%.0f kip-ft) at x = %.0f in\n', ...
    results.M(idx_M_max), results.M(idx_M_max)/12, results.x(idx_M_max));
fprintf('  V_max = %.1f kips at x = %.0f in\n', ...
    results.V(idx_V_max), results.x(idx_V_max));
fprintf('  N_max = %.1f kips (prestress compression)\n', N_max);

fprintf('\nPrestress Effects at Midspan:\n');
mid_idx = round(length(results.x)/2);
fprintf('  Effective prestress force: P = %.1f kips\n', results.P(mid_idx));
fprintf('  Eccentricity: e = %.2f in\n', results.e(mid_idx));
fprintf('  Prestress moment: M_p = %.0f kip-in\n', results.M_prestress(mid_idx));

fprintf('\nStresses at Midspan:\n');
fprintf('  Top fiber (total): f_t = %.3f ksi\n', results.stresses.f_top_total(mid_idx));
fprintf('  Bottom fiber (total): f_b = %.3f ksi\n', results.stresses.f_bot_total(mid_idx));
fprintf('  Allowable compression: %.3f ksi\n', results.stresses.fc_allow_compression);
fprintf('  Allowable tension: %.3f ksi\n', results.stresses.fc_allow_tension);

if isfield(results, 'capacity')
    fprintf('\nNominal Capacity:\n');
    fprintf('  M_n = %.0f kip-in (%.0f kip-ft)\n', results.capacity.Mn, results.capacity.Mn/12);
    fprintf('  phi*M_n = %.0f kip-in (%.0f kip-ft)\n', results.capacity.phi_Mn, results.capacity.phi_Mn/12);
    fprintf('  Demand/Capacity Ratio: %.2f\n', abs(results.M(idx_M_max))/results.capacity.phi_Mn);
end

fprintf('\nReactions:\n');
if isfield(results.reactions, 'left')
    fprintf('  R_left = %.1f kips at x = %.0f in\n', results.reactions.left, results.reactions.x_left);
    fprintf('  R_right = %.1f kips at x = %.0f in\n', results.reactions.right, results.reactions.x_right);
end

%% Generate plots
fprintf('\nGenerating plots...\n');

% Main results plot (N, V, M diagrams)
options.show_stresses = true;
options.show_section = true;
plotPrestressedBeamResults(results, options);

% Cross-section at midspan
section_options.beam = beam;
section_options.x_location = 0.5;  % Midspan
section_options.show_dimensions = true;
section_options.show_properties = true;
plotSection(section, prestress, reinforcement, materials, section_options);

% Cross-section at support (different tendon position)
figure;
section_options.x_location = 0.0;  % At support
section_options.show_dimensions = true;
plotSection(section, prestress, reinforcement, materials, section_options);
title('Cross-Section at Support (x = 0)', 'FontWeight', 'bold');


% At midspan (using fraction)
plotSectionStressStrain(results, 0.5);

% At specific location (in inches)
plotSectionStressStrain(results, 300);

% With options
options.plot_type = 'both';        % 'stress', 'strain', or 'both'
options.show_components = true;    % Show prestress/external separately
options.show_limits = true;        % Show allowable stress limits
plotSectionStressStrain(results, 0.5, options);

%% Export results to workspace
fprintf('\nResults saved to workspace variable: results\n');

%% Optional: Generate detailed report
generateReport = false;  % Set to true to generate text report
if generateReport
    generateAnalysisReport(results, 'PrestressedBeamReport.txt');
end

fprintf('\n========================================\n');
fprintf('Analysis complete!\n');
fprintf('========================================\n');
