%% MAIN_PRESTRESSED_BEAM_ANALYSIS
% Main driver script for prestressed concrete beam analysis
% Runs multiple ACI 318-19 design stages (Transfer, Service-Sustained, Service-Total)
% and saves figures per stage to the output folder.
%
% Usage:
%   1. Edit inputPrestressedBeam_Project1.m to define your beam
%   2. Edit defineDesignStages.m to toggle which stages to compute
%   3. Run this script — figures are saved to Project1/<StageName>_Figure_*.png

clear; clc; close all;
My_data = 'Project1';

%% ---- Load input data (single input file for all stages) ----
fprintf('Loading input data...\n');
[beam, section, materials, prestress, reinforcement, loads] = inputPrestressedBeam_Project1();

%% ---- Define design stages ----
%   Edit defineDesignStages.m to toggle stages on/off (stage.active = false to skip)
stages = defineDesignStages(materials, prestress);

%% ---- Cross-section plot locations — defined in inputPrestressedBeam_Project1.m ----
x_fractions = beam.x_plot_fractions;

%% ========================================================================
%%  STAGE LOOP
%% ========================================================================
results_all = {};

for i = 1:numel(stages)
    stg = stages{i};

    if ~stg.active
        fprintf('\n[SKIP] Stage: %s  (active = false)\n', stg.name);
        continue;
    end

    fprintf('\n========================================\n');
    fprintf('  STAGE: %s\n', stg.name);
    fprintf('  eta = %.2f  |  SW=%d  SDL=%d  LL=%d\n', ...
        stg.eta, stg.include_sw, stg.include_sdl, stg.include_ll);
    fprintf('  Allowable:  compr = +%.3f ksi  |  tens = %.3f ksi\n', ...
        stg.f_allow_compr, stg.f_allow_tens);
    fprintf('========================================\n');

    %% Run analysis for this stage
    results_i = analyzePrestressedBeam(beam, section, materials, prestress, ...
                                        reinforcement, loads, stg);
    results_all{end+1} = results_i;

    %% Print summary
    printStageSummary(results_i);

    %% Generate plots
    close all;   % reset figure numbers for each stage

    % --- Main N/V/M diagram ---
    options.show_stresses = true;
    options.show_section  = true;
    options.stage_name    = stg.name;
    plotPrestressedBeamResults(results_i, options);

    % --- Cross-section geometry at each location ---
    for k = 1:length(x_fractions)
        sec_opts.beam             = beam;
        sec_opts.x_location       = x_fractions(k);
        sec_opts.show_dimensions  = true;
        sec_opts.show_properties  = true;
        plotSection(section, prestress, reinforcement, materials, sec_opts);
    end

    % --- Stress/strain distribution at each location ---
    for k = 1:length(x_fractions)
        ss_opts.stage_label = stg.name;
        plotSectionStressStrain(results_i, x_fractions(k), ss_opts);
    end

    %% Save all figures for this stage
    saveStageFigures(My_data, stg.name);
end

fprintf('\n========================================\n');
fprintf('All active stages complete.\n');
fprintf('Figures saved to: %s\\\n', My_data);
fprintf('========================================\n');

%% Export last results to workspace
if ~isempty(results_all)
    results = results_all{end};   %#ok<NASGU>
    fprintf('Last stage results saved to workspace variable: results\n');
end


%% ========================================================================
%%  LOCAL FUNCTIONS
%% ========================================================================

function printStageSummary(results_i)
% Print key results for a single stage to the command window.

fprintf('\n--- Results Summary: %s ---\n', results_i.stage_name);

[M_max, idx_M_max] = max(abs(results_i.M));
[V_max, idx_V_max] = max(abs(results_i.V));
[N_max, ~]         = max(abs(results_i.N));
mid_idx            = round(length(results_i.x) / 2);

fprintf('\n  Maximum Internal Forces:\n');
fprintf('    M_max = %.0f kip-in (%.0f kip-ft) at x = %.0f in\n', ...
    results_i.M(idx_M_max), results_i.M(idx_M_max)/12, results_i.x(idx_M_max));
fprintf('    V_max = %.1f kips at x = %.0f in\n', ...
    results_i.V(idx_V_max), results_i.x(idx_V_max));
fprintf('    N_max = %.1f kips (prestress compression)\n', N_max);

fprintf('\n  Prestress Effects at Midspan:\n');
fprintf('    Effective prestress force: P = %.1f kips\n', results_i.P(mid_idx));
fprintf('    Eccentricity: e = %.2f in\n', results_i.e(mid_idx));
fprintf('    Prestress moment: M_p = %.0f kip-in\n', results_i.M_prestress(mid_idx));

fprintf('\n  Stresses at Midspan:\n');
fprintf('    Top fiber:  f_t = %+.3f ksi  (allow: [%.3f, +%.3f])\n', ...
    results_i.stresses.f_top_total(mid_idx), ...
    results_i.stresses.fc_allow_tension, results_i.stresses.fc_allow_compression);
fprintf('    Bot fiber:  f_b = %+.3f ksi  (allow: [%.3f, +%.3f])\n', ...
    results_i.stresses.f_bot_total(mid_idx), ...
    results_i.stresses.fc_allow_tension, results_i.stresses.fc_allow_compression);

if isfield(results_i, 'capacity')
    fprintf('\n  Nominal Capacity:\n');
    fprintf('    M_n     = %.0f kip-in (%.0f kip-ft)\n', ...
        results_i.capacity.Mn, results_i.capacity.Mn/12);
    fprintf('    φ·M_n   = %.0f kip-in (%.0f kip-ft)\n', ...
        results_i.capacity.phi_Mn, results_i.capacity.phi_Mn/12);
    fprintf('    D/C     = %.2f\n', abs(results_i.M(idx_M_max)) / results_i.capacity.phi_Mn);
end

if isfield(results_i.reactions, 'left')
    fprintf('\n  Reactions:\n');
    fprintf('    R_left  = %.1f kips at x = %.0f in\n', ...
        results_i.reactions.left,  results_i.reactions.x_left);
    fprintf('    R_right = %.1f kips at x = %.0f in\n', ...
        results_i.reactions.right, results_i.reactions.x_right);
end
end


function saveStageFigures(output_dir, stage_name)
% Save all open figures with stage-prefixed filenames to output_dir.

fig_handles = findall(0, 'Type', 'figure');
if isempty(fig_handles);  return;  end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = 1:length(fig_handles)
    fig_num  = fig_handles(i).Number;
    fig_name = sprintf('%s_Figure_%d', stage_name, fig_num);

    % Try to append descriptive name from figure title
    if ~isempty(fig_handles(i).Children)
        try
            title_obj = get(fig_handles(i).Children(1), 'Title');
            if ~isempty(title_obj) && ~isempty(title_obj.String)
                title_str = title_obj.String;
                if iscell(title_str);  title_str = title_str{1};  end
                title_str = regexprep(title_str, '[^\w\s-]', '');
                title_str = strtrim(regexprep(title_str, '\s+', '_'));
                if ~isempty(title_str)
                    fig_name = sprintf('%s_Figure_%d_%s', stage_name, fig_num, title_str);
                end
            end
        catch
            % keep default name
        end
    end

    figure(fig_handles(i));
    saveas(fig_handles(i), fullfile(output_dir, [fig_name '.png']));
    savefig(fig_handles(i), fullfile(output_dir, [fig_name '.fig']));
    fprintf('  Saved: %s\n', fig_name);
end

fprintf('  → %d figure(s) saved to %s\\\n', length(fig_handles), output_dir);
end
