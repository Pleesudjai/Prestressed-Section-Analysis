%% RUN_END_BLOCK_DESIGN.m
% End Block (Anchorage Zone) Design — Gergely-Sozen Method
% Standalone script for the Double-T prestressed beam.
%
% HOW TO RUN:
%   Set MATLAB current folder to 08_matlab_program_analysis/, then >> runEndBlockDesign
%
% OUTPUTS:
%   - Console: full calculation detail (stresses, net moments, stirrup sizing)
%   - Figures: cross-section, stress/strain at x=0, net moment diagram,
%              reinforcement layout — saved to project_Project2/output/EndBlock/
%
% METHOD:
%   Gergely-Sozen free-body approach (Naaman Sec. 4.17)
%   T_burst = M_max / (3h/4), f_s = 20 ksi (WSD)
%
% SIGN CONVENTION: compression = positive, tension = negative
% COORDINATE ORIGIN: y = 0 at bottom of stems

clear; clc; close all;

output_dir = fullfile('project_Project2', 'output', 'EndBlock');

%% ================================================================
%%  1. LOAD INPUT DATA
%% ================================================================
fprintf('Loading input data...\n');
[beam, section, materials, prestress, reinforcement, loads] = inputData();

%% ================================================================
%%  2. RUN END BLOCK DESIGN
%% ================================================================
eb = endBlockDesign(beam, section, materials, prestress);

%% ================================================================
%%  3. CROSS-SECTION AND STRESS PLOTS (reuse existing templates)
%% ================================================================
% Run Transfer stage to get results at x = 0
stages = defineDesignStages(materials, prestress);
for i = 1:numel(stages)
    if contains(stages{i}.name, 'Transfer', 'IgnoreCase', true) && stages{i}.active
        results_transfer = analyzePrestressedBeam(beam, section, materials, ...
            prestress, reinforcement, loads, stages{i});
        break;
    end
end

if exist('results_transfer', 'var')
    % Cross-section at support
    sec_opts.beam             = beam;
    sec_opts.x_location       = 0;
    sec_opts.show_dimensions  = true;
    sec_opts.show_properties  = true;
    plotSection(section, prestress, reinforcement, materials, sec_opts);

    % Stress/strain distribution at x = 0
    ss_opts.stage_label = 'EndBlock - Transfer';
    plotSectionStressStrain(results_transfer, 0, ss_opts);
end

%% ================================================================
%%  4. SAVE ALL FIGURES
%% ================================================================
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

fig_handles = findall(0, 'Type', 'figure');
for i = 1:length(fig_handles)
    if ~isvalid(fig_handles(i)); continue; end
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
    try
        savefig(fig_handles(i), fullfile(output_dir, [name_str '.fig']));
    catch
    end
    fprintf('Saved: %s\n', name_str);
end
fprintf('Total: %d figures saved to %s/\n', length(fig_handles), output_dir);

%% Export to workspace
fprintf('\nEnd block results saved to workspace variable: eb\n');
