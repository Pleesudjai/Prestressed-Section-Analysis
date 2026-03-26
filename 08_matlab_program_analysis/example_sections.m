%% EXAMPLE_SECTIONS - Demonstrate various section geometries
% Shows how to define different prestressed concrete sections
%
% Run this script to see examples of section definitions

clear; clc; close all;

fprintf('===========================================\n');
fprintf('  PRESTRESSED BEAM SECTION EXAMPLES\n');
fprintf('===========================================\n\n');

%% Example 1: Simple Rectangle
fprintf('Example 1: Rectangular Section\n');
fprintf('------------------------------\n');

rect_vertices = [
    0,  0;      % Bottom-left
    18, 0;      % Bottom-right
    18, 36;     % Top-right
    0,  36;     % Top-left
];

section_rect = createPolygonalSection(rect_vertices, 'Name', 'Rectangle 18x36');

%% Example 2: I-Section (AASHTO Type)
fprintf('Example 2: I-Section\n');
fprintf('--------------------\n');

I_vertices = [
    0,    0;        % Bottom flange
    24,   0;
    24,   6;
    15,   6;        % Web transition
    15,   30;
    24,   30;       % Top flange
    24,   36;
    0,    36;
    0,    30;
    9,    30;       % Web transition
    9,    6;
    0,    6;
];

section_I = createPolygonalSection(I_vertices, 'Name', 'I-Section');

%% Example 3: T-Section
fprintf('Example 3: T-Section\n');
fprintf('--------------------\n');

T_vertices = [
    0,    0;        % Web bottom
    12,   0;
    12,   24;       % Flange bottom
    0,    24;       % Web-flange junction
    0,    24;       % Return to web (need to add flange)
    -12,  24;       % Flange extends left
    -12,  30;       % Flange top
    24,   30;       % Flange extends right
    24,   24;
    12,   24;
    12,   0;        % Back to start
];

% Corrected T-section
T_vertices = [
    6,    0;        % Web left-bottom
    18,   0;        % Web right-bottom
    18,   24;       % Web right-top
    36,   24;       % Flange right-bottom
    36,   30;       % Flange right-top
    0,    30;       % Flange left-top
    0,    24;       % Flange left-bottom
    6,    24;       % Web left-top
];

section_T = createPolygonalSection(T_vertices, 'Name', 'T-Section');

%% Example 4: Box Section (with hole)
fprintf('Example 4: Box Section (Hollow)\n');
fprintf('-------------------------------\n');

box_outer = [
    0,    0;
    36,   0;
    36,   30;
    0,    30;
];

box_inner = [
    6,    6;        % Clockwise for hole
    6,    24;
    30,   24;
    30,   6;
];

section_box = createPolygonalSection(box_outer, 'Holes', {box_inner}, 'Name', 'Box Section');

%% Example 5: L-Section (similar to inputData.m)
fprintf('Example 5: L-Section\n');
fprintf('--------------------\n');

L_vertices = [
    20,   240;      % Top of vertical leg
    0,    240;
    0,    0;        % Origin
    150,  0;        % Horizontal leg extends right
    150,  20;
    20,   20;       % Inside corner
];

section_L = createPolygonalSection(L_vertices, 'Name', 'L-Section');

%% Example 6: Double-T Section
fprintf('Example 6: Double-T Section\n');
fprintf('---------------------------\n');

DT_vertices = [
    0,    0;        % Left stem bottom
    4,    0;
    4,    18;       % Left stem top
    20,   18;       % Between stems
    20,   0;        % Right stem bottom
    24,   0;
    24,   18;       % Right stem top
    36,   18;       % Flange overhang right
    36,   24;       % Flange top
    0,    24;       % Flange top left
    0,    18;       % Flange bottom left
    0,    0;        % Close
];

% Corrected Double-T
DT_vertices = [
    0,    0;        % Left stem left
    6,    0;        % Left stem right
    6,    18;       % Left stem top right
    18,   18;       % Between stems bottom
    18,   0;        % Right stem left
    24,   0;        % Right stem right
    24,   18;       % Right stem top left
    36,   18;       % Flange right bottom
    36,   24;       % Flange right top
    -12,  24;       % Flange left top
    -12,  18;       % Flange left bottom
    0,    18;       % Back to left stem
];

section_DT = createPolygonalSection(DT_vertices, 'Name', 'Double-T Section');

%% Example 7: Trapezoidal Section
fprintf('Example 7: Trapezoidal Section\n');
fprintf('------------------------------\n');

trap_vertices = [
    6,    0;        % Bottom-left
    30,   0;        % Bottom-right
    24,   24;       % Top-right
    12,   24;       % Top-left
];

section_trap = createPolygonalSection(trap_vertices, 'Name', 'Trapezoidal Section');

%% Example 8: Bulb-T Section (Simplified)
fprintf('Example 8: Bulb-T Section (Simplified)\n');
fprintf('--------------------------------------\n');

bulbT_vertices = [
    6,    0;        % Bulb left
    18,   0;        % Bulb right
    16,   8;        % Bulb-web transition right
    14,   8;        % Web bottom right
    14,   30;       % Web top right
    36,   30;       % Flange right
    36,   36;       % Flange top right
    0,    36;       % Flange top left
    0,    30;       % Flange left
    10,   30;       % Web top left
    10,   8;        % Web bottom left
    8,    8;        % Bulb-web transition left
];

section_bulbT = createPolygonalSection(bulbT_vertices, 'Name', 'Bulb-T Section');

%% Visualize all sections
fprintf('\n===========================================\n');
fprintf('  VISUALIZING ALL SECTIONS\n');
fprintf('===========================================\n');

figure('Name', 'Section Gallery', 'Position', [50, 50, 1400, 800], 'Color', 'w');

sections = {section_rect, section_I, section_T, section_box, section_L, section_DT, section_trap, section_bulbT};
titles = {'Rectangle', 'I-Section', 'T-Section', 'Box Section', 'L-Section', 'Double-T', 'Trapezoidal', 'Bulb-T'};

for i = 1:length(sections)
    subplot(2, 4, i);
    hold on;
    axis equal;
    grid on;
    box on;
    
    sect = sections{i};
    v = sect.vertices;
    v_closed = [v; v(1,:)];
    
    % Fill section
    fill(v(:,1), v(:,2), [0.7, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5);
    
    % Plot holes if present
    if isfield(sect, 'holes') && ~isempty(sect.holes)
        for j = 1:length(sect.holes)
            hole = sect.holes{j};
            fill(hole(:,1), hole(:,2), 'w', 'EdgeColor', 'k', 'LineWidth', 1);
        end
    end
    
    % Mark centroid
    plot(sect.xc, sect.yc, 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    
    title(sprintf('%s\nA=%.0f, I_x=%.0f', titles{i}, sect.A, sect.Ix), 'FontSize', 10);
    xlabel('x (in)');
    ylabel('y (in)');
end

sgtitle('Prestressed Concrete Section Gallery', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nSection gallery displayed.\n');

%% Quick Analysis Example
fprintf('\n===========================================\n');
fprintf('  QUICK ANALYSIS EXAMPLE\n');
fprintf('===========================================\n');

% Use I-section with P = 300 kips, parabolic e profile
P = 300;  % kips
e = [4; 10; 4];  % eccentricity at start, mid, end (in)
L = 480;  % 40 ft span

fprintf('\nRunning quick analysis with:\n');
fprintf('  P = %.0f kips\n', P);
fprintf('  e = [%.1f, %.1f, %.1f] in (parabolic)\n', e(1), e(2), e(3));
fprintf('  L = %.0f in (%.1f ft)\n', L, L/12);

[beam, section, materials, prestress, reinforcement, loads] = ...
    quickInputPE(P, e, L, I_vertices, ...
    'ProfileType', 'parabolic', ...
    'DistributedLoad', -0.05, ...
    'PointLoads', [L/3, -20, 0; 2*L/3, -20, 0]);

% Run analysis
results = analyzePrestressedBeam(beam, section, materials, prestress, reinforcement, loads);

% Plot results
plotPrestressedBeamResults(results);

fprintf('\nAnalysis complete! Check figures for results.\n');
