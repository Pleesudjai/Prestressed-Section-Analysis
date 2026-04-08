%% MAIN_PRESTRESSED_BEAM_ANALYSIS
% Main driver script — runs inputData.m in the current project folder.
% Shared helper functions are loaded from the parent directory.
%
% HOW TO RUN:
%   1. Set MATLAB current folder to this project folder (project_xxx/)
%   2. Run this script  (F5 or >> main_PrestressedBeamAnalysis)
%   3. Figures are saved to  output/<StageName>_Figure_*.png
%
% DESIGN STAGES:
%   Edit defineDesignStages.m (in parent folder) to toggle stages on/off.
%   Default: Transfer  |  Service_Sustained  |  Service_Total

%% Add parent folder so shared functions are found
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

clear; clc; close all;
My_data = 'output';   % figures saved to  <project_folder>/output/

%% Load input data
fprintf('Loading input data...\n');
[beam, section, materials, prestress, reinforcement, loads] = inputData();

%% Print allowable stresses
fprintf('\n----------------------------------------\n');
fprintf('  ALLOWABLE STRESSES  [Edition: %s]\n', materials.code_edition);
fprintf('----------------------------------------\n');
fprintf('  AT TRANSFER (f''ci = %.1f ksi):\n', materials.fci);
fprintf('    Compression (general) : +%.3f ksi  (+0.60 f''ci)\n', materials.f_ci_allow);
fprintf('    Compression (ends)    : +%.3f ksi  (+%.2f f''ci)\n', ...
    materials.f_ci_allow_end, materials.f_ci_allow_end / materials.fci);
fprintf('    Tension (general)     : %.3f ksi  (-3*sqrt(f''ci))\n', materials.f_ti_allow);
fprintf('    Tension (ends)        : %.3f ksi  (-6*sqrt(f''ci))\n', materials.f_ti_allow_end);
fprintf('  AT SERVICE (f''c = %.1f ksi):\n', materials.fc);
fprintf('    Compression (sustained SW+SDL) : +%.3f ksi  (+0.45 f''c)\n', materials.f_cs_allow_sust);
fprintf('    Compression (total SW+SDL+LL)  : +%.3f ksi  (+0.60 f''c)\n', materials.f_cs_allow_total);
fprintf('    Tension Class U boundary       : %.3f ksi  (-%.1f*sqrt(f''c))\n', ...
    materials.f_tu_allow, abs(materials.f_tu_allow) / sqrt(materials.fc * 1000) * 1000);
fprintf('    Tension (Class C)              : %.3f ksi  (-12*sqrt(f''c))\n', materials.f_ts_allow);
fprintf('----------------------------------------\n');

%% Define design stages
%   Edit defineDesignStages.m in the parent folder to toggle stages on/off
stages = defineDesignStages(materials, prestress);

%% Cross-section plot locations — defined in inputData.m (beam.x_plot_fractions)
x_fractions = beam.x_plot_fractions;   % fractions of span, e.g. [0, 0.25, 0.50, 0.75]

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

    % Main N/V/M diagram
    options.show_stresses = true;
    options.show_section  = true;
    options.stage_name    = stg.name;
    plotPrestressedBeamResults(results_i, options);

    % Cross-section geometry at each location
    for k = 1:length(x_fractions)
        sec_opts.beam             = beam;
        sec_opts.x_location       = x_fractions(k);
        sec_opts.show_dimensions  = true;
        sec_opts.show_properties  = true;
        plotSection(section, prestress, reinforcement, materials, sec_opts);
    end

    % Stress/strain distribution at each location
    for k = 1:length(x_fractions)
        ss_opts.stage_label = stg.name;
        plotSectionStressStrain(results_i, x_fractions(k), ss_opts);
    end

    %% Save all figures for this stage
    saveStageFigures(My_data, stg.name);
end

fprintf('\n========================================\n');
fprintf('All active stages complete.\n');
fprintf('Figures saved to: %s/\n', My_data);
fprintf('========================================\n');

%% Export last results to workspace
if ~isempty(results_all)
    results = results_all{end};   %#ok<NASGU>
    fprintf('Last stage results saved to workspace variable: results\n');
end

%% End block design (Gergely-Sozen method)
fprintf('\n========================================\n');
fprintf('Running end block design...\n');
eb = endBlockDesign(beam, section, materials, prestress);
saveStageFigures(My_data, 'EndBlock');
fprintf('========================================\n');

%% Generate calculation report (.txt)
fprintf('\n========================================\n');
fprintf('Generating calculation report...\n');
generateReport(beam, section, materials, prestress, loads, results_all, My_data);
fprintf('========================================\n');


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
fprintf('    P = %.1f kips | e = %.2f in | M_p = %.0f kip-in\n', ...
    results_i.P(mid_idx), results_i.e(mid_idx), results_i.M_prestress(mid_idx));

fprintf('\n  Stresses at Midspan:\n');
fprintf('    f_top = %+.3f ksi  (allow: [%.3f, +%.3f])\n', ...
    results_i.stresses.f_top_total(mid_idx), ...
    results_i.stresses.fc_allow_tension, results_i.stresses.fc_allow_compression);
fprintf('    f_bot = %+.3f ksi  (allow: [%.3f, +%.3f])\n', ...
    results_i.stresses.f_bot_total(mid_idx), ...
    results_i.stresses.fc_allow_tension, results_i.stresses.fc_allow_compression);

if isfield(results_i, 'capacity')
    fprintf('\n  Nominal Capacity:\n');
    fprintf('    M_n   = %.0f kip-in (%.0f kip-ft)\n', ...
        results_i.capacity.Mn, results_i.capacity.Mn/12);
    fprintf('    phi*Mn = %.0f kip-in (%.0f kip-ft)\n', ...
        results_i.capacity.phi_Mn, results_i.capacity.phi_Mn/12);
    fprintf('    D/C   = %.2f\n', abs(results_i.M(idx_M_max)) / results_i.capacity.phi_Mn);
end

if isfield(results_i.reactions, 'left')
    fprintf('\n  Reactions:\n');
    fprintf('    R_left  = %.1f kips | R_right = %.1f kips\n', ...
        results_i.reactions.left, results_i.reactions.right);
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
        end
    end

    if ~isvalid(fig_handles(i));  continue;  end
    try
        figure(fig_handles(i));
        saveas(fig_handles(i), fullfile(output_dir, [fig_name '.png']));
        try
            savefig(fig_handles(i), fullfile(output_dir, [fig_name '.fig']));
        catch
            % savefig can fail in batch/headless mode — skip .fig silently
        end
        fprintf('  Saved: %s\n', fig_name);
    catch ME
        fprintf('  Warning: could not save %s: %s\n', fig_name, ME.message);
    end
end

fprintf('  → %d figure(s) saved to %s/\n', length(fig_handles), output_dir);
end
