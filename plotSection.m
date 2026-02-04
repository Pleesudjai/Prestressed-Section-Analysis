function plotSection(section, prestress, reinforcement, materials, options)
% PLOTSECTION - Visualize prestressed concrete cross-section with detailed dimensions
%
% Inputs:
%   section - Section geometry structure
%   prestress - Prestressing tendon information
%   reinforcement - Non-prestressed reinforcement
%   materials - Material properties
%   options - Optional visualization options
%
% Usage:
%   plotSection(section, prestress, reinforcement, materials)
%
% Features:
%   - No grid background for clean drawings
%   - Comprehensive dimensions for all flanges and web
%   - Rebar and tendon location dimensions from bottom
%   - Suitable for copying section drawings

if nargin < 5
    options = struct();
end

% Default options
if ~isfield(options, 'show_dimensions'), options.show_dimensions = true; end
if ~isfield(options, 'show_properties'), options.show_properties = false; end
if ~isfield(options, 'x_location'), options.x_location = 0.5; end
if ~isfield(options, 'beam'), options.beam = []; end
if ~isfield(options, 'show_centroid'), options.show_centroid = true; end

%% Create figure
figure('Name', 'Cross-Section', 'Position', [100, 100, 900, 800], 'Color', 'w');

hold on;
axis equal;
grid off;  % NO GRID
box off;   % Clean look

%% Plot concrete section
v = section.vertices;

% Fill concrete with light gray
fill(v(:,1), v(:,2), [0.85, 0.85, 0.85], 'FaceAlpha', 0.7, 'EdgeColor', 'k', 'LineWidth', 2);

%% Extract section geometry for dimensioning
[sectionType, dims] = analyzeSectionGeometry(v);

%% Plot centroid and centroidal axes
xc = section.xc;
yc = section.yc;
y_bottom = min(v(:,2));  % Bottom of section for measuring rebar locations

if options.show_centroid
    plot(xc, yc, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
    
    % Centroidal axes (dashed)
    x_range = max(v(:,1)) - min(v(:,1));
    y_range = max(v(:,2)) - min(v(:,2));
    axis_ext = max(x_range, y_range) * 0.15;
    
    % Horizontal centroidal axis
    plot([min(v(:,1)) - axis_ext*0.5, max(v(:,1)) + axis_ext*0.5], [yc, yc], ...
        'k--', 'LineWidth', 0.8);
    % Vertical centroidal axis  
    plot([xc, xc], [min(v(:,2)) - axis_ext*0.3, max(v(:,2)) + axis_ext*0.3], ...
        'k--', 'LineWidth', 0.8);
    
    % Label centroid
    text(xc + 2, yc + 2, 'C.G.', 'FontSize', 10, 'FontWeight', 'bold');
end

%% Collect rebar/tendon y-positions for dimensioning
rebar_positions = [];  % [y_position, label, x_position]

%% Plot non-prestressed reinforcement
if ~isempty(reinforcement) && isfield(reinforcement, 'longitudinal')
    for i = 1:length(reinforcement.longitudinal)
        rebar = reinforcement.longitudinal{i};
        if isempty(rebar.y_position) || isempty(rebar.x_positions)
            continue;
        end
        y_pos = rebar.y_position;
        x_pos = rebar.x_positions;
        areas = rebar.areas;
        
        % Store position for dimensioning
        rebar_positions = [rebar_positions; y_pos, mean(x_pos), sum(areas)];
        
        for j = 1:length(x_pos)
            d_bar = 2 * sqrt(areas(j) / pi);
            theta = linspace(0, 2*pi, 30);
            x_circle = x_pos(j) + (d_bar/2) * cos(theta);
            y_circle = y_pos + (d_bar/2) * sin(theta);
            fill(x_circle, y_circle, 'k');
        end
    end
end

%% Plot prestressing tendons
tendon_positions = [];  % [y_position, x_position, Aps]

if ~isempty(prestress) && isfield(prestress, 'tendons')
    colors = lines(length(prestress.tendons));
    
    if ~isempty(options.beam)
        L = options.beam.L;
        x_loc = options.x_location * L;
        
        for i = 1:length(prestress.tendons)
            tendon = prestress.tendons{i};
            e = interp1(options.beam.x, tendon.e, x_loc, 'linear');
            y_tendon = yc - e;
            d_tendon = 2 * sqrt(tendon.Aps / pi);
            bonded = interp1(options.beam.x, double(tendon.bonded), x_loc, 'nearest');
            
            % Store position for dimensioning
            tendon_positions = [tendon_positions; y_tendon, xc, tendon.Aps, i];
            
            theta = linspace(0, 2*pi, 30);
            x_circle = xc + (d_tendon/2) * cos(theta);
            y_circle = y_tendon + (d_tendon/2) * sin(theta);
            
            if bonded
                fill(x_circle, y_circle, colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1);
            else
                fill(x_circle, y_circle, 'w', 'EdgeColor', colors(i,:), 'LineWidth', 2);
            end
        end
    else
        for i = 1:length(prestress.tendons)
            tendon = prestress.tendons{i};
            e = tendon.e_mid;
            y_tendon = yc - e;
            d_tendon = 2 * sqrt(tendon.Aps / pi);
            
            % Store position for dimensioning
            tendon_positions = [tendon_positions; y_tendon, xc, tendon.Aps, i];
            
            theta = linspace(0, 2*pi, 30);
            x_circle = xc + (d_tendon/2) * cos(theta);
            y_circle = y_tendon + (d_tendon/2) * sin(theta);
            
            fill(x_circle, y_circle, colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1);
        end
    end
end

%% Add comprehensive dimensions
if options.show_dimensions
    drawComprehensiveDimensions(v, section, sectionType, dims, rebar_positions, tendon_positions, y_bottom);
end

%% Axis labels and title
xlabel('x (in)', 'FontSize', 11);
ylabel('y (in)', 'FontSize', 11);

if ~isempty(options.beam) && options.x_location > 0
    title(sprintf('Cross-Section at x = %.0f in (%.1f%% of span)', ...
        options.x_location * options.beam.L, options.x_location * 100), ...
        'FontWeight', 'bold', 'FontSize', 12);
else
    title('Cross-Section', 'FontWeight', 'bold', 'FontSize', 12);
end

% NO LEGEND - removed as requested

% Adjust axis limits to show all dimensions
axis auto;
curr_xlim = xlim();
curr_ylim = ylim();
xlim([curr_xlim(1) - 5, curr_xlim(2) + 5]);
ylim([curr_ylim(1) - 5, curr_ylim(2) + 5]);

set(gca, 'FontSize', 10);

end

%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function [sectionType, dims] = analyzeSectionGeometry(v)
% Analyze polygon vertices to determine section type and extract key dimensions

dims = struct();
n = size(v, 1);

% Basic bounds
dims.x_min = min(v(:,1));
dims.x_max = max(v(:,1));
dims.y_min = min(v(:,2));
dims.y_max = max(v(:,2));
dims.total_width = dims.x_max - dims.x_min;
dims.total_height = dims.y_max - dims.y_min;

% Find unique y-coordinates (horizontal transitions)
y_unique = unique(v(:,2));
y_unique = sort(y_unique);

% Find unique x-coordinates (vertical transitions)  
x_unique = unique(v(:,1));
x_unique = sort(x_unique);

% Determine section type based on number of y-levels
if length(y_unique) == 4
    % Likely I-section or T-section
    x_center = (dims.x_min + dims.x_max) / 2;
    
    % Get widths at each y-level
    widths_at_y = zeros(length(y_unique), 2);  % [x_left, x_right]
    for i = 1:length(y_unique)
        y_level = y_unique(i);
        tol = 0.01;
        idx = abs(v(:,2) - y_level) < tol;
        x_at_level = v(idx, 1);
        if ~isempty(x_at_level)
            widths_at_y(i,:) = [min(x_at_level), max(x_at_level)];
        end
    end
    
    sectionType = 'I-section';
    
    % Bottom and top extremes
    dims.bot_flange_width = widths_at_y(1,2) - widths_at_y(1,1);
    dims.bot_flange_height = y_unique(2) - y_unique(1);
    dims.bot_flange_y = [y_unique(1), y_unique(2)];
    dims.bot_flange_x = widths_at_y(1,:);
    
    dims.top_flange_width = widths_at_y(4,2) - widths_at_y(4,1);
    dims.top_flange_height = y_unique(4) - y_unique(3);
    dims.top_flange_y = [y_unique(3), y_unique(4)];
    dims.top_flange_x = widths_at_y(4,:);
    
    dims.web_height = y_unique(3) - y_unique(2);
    dims.web_y = [y_unique(2), y_unique(3)];
    
    % Find web width by looking at the INNER x-coordinates at transition levels
    y_bot_transition = y_unique(2);
    y_top_transition = y_unique(3);
    
    idx_bot = abs(v(:,2) - y_bot_transition) < 0.01;
    idx_top = abs(v(:,2) - y_top_transition) < 0.01;
    
    x_at_bot = sort(v(idx_bot, 1));
    x_at_top = sort(v(idx_top, 1));
    
    if length(x_at_bot) >= 4
        web_x_left = x_at_bot(2);
        web_x_right = x_at_bot(end-1);
    elseif length(x_at_bot) == 2
        web_x_left = x_at_bot(1);
        web_x_right = x_at_bot(2);
    else
        web_x_left = min(x_at_bot);
        web_x_right = max(x_at_bot);
    end
    
    dims.web_width = web_x_right - web_x_left;
    dims.web_x = [web_x_left, web_x_right];
    
elseif length(y_unique) == 2
    % Rectangular section
    sectionType = 'rectangular';
    dims.width = dims.total_width;
    dims.height = dims.total_height;
    
else
    % General polygon (L-section, T-section, etc.)
    sectionType = 'general';
end

dims.y_levels = y_unique;
dims.x_levels = x_unique;

end

function drawComprehensiveDimensions(v, section, sectionType, dims, rebar_positions, tendon_positions, y_bottom)
% Draw comprehensive dimensions for the section

% Dimension line parameters - scale with section size
scale = max(dims.total_width, dims.total_height);
offset = scale * 0.12;
ext = scale * 0.03;
gap = scale * 0.02;
fontsize = 9;

switch sectionType
    case 'I-section'
        drawIBeamDimensions(v, section, dims, offset, ext, gap, fontsize);
    case 'rectangular'
        drawRectangularDimensions(v, section, dims, offset, ext, gap, fontsize);
    otherwise
        drawGeneralDimensions(v, section, dims, offset, ext, gap, fontsize);
end

% Draw rebar and tendon location dimensions
drawRebarTendonDimensions(dims, rebar_positions, tendon_positions, y_bottom, offset, ext, fontsize);

end

function drawRebarTendonDimensions(dims, rebar_positions, tendon_positions, y_bottom, offset, ext, fontsize)
% Draw dimension lines showing rebar and tendon locations from bottom

% Use a position further left for rebar/tendon dimensions
x_dim_rebar = dims.x_min - offset * 1.8;

% Combine all positions and sort by y
all_positions = [];

% Add rebar positions
if ~isempty(rebar_positions)
    for i = 1:size(rebar_positions, 1)
        y_pos = rebar_positions(i, 1);
        As = rebar_positions(i, 3);
        all_positions = [all_positions; y_pos, 1, As, 0];  % type=1 for rebar
    end
end

% Add tendon positions
if ~isempty(tendon_positions)
    for i = 1:size(tendon_positions, 1)
        y_pos = tendon_positions(i, 1);
        Aps = tendon_positions(i, 3);
        tendon_num = tendon_positions(i, 4);
        all_positions = [all_positions; y_pos, 2, Aps, tendon_num];  % type=2 for tendon
    end
end

if isempty(all_positions)
    return;
end

% Sort by y-position
[~, sort_idx] = sort(all_positions(:,1));
all_positions = all_positions(sort_idx, :);

% Draw dimension from bottom to each rebar/tendon
for i = 1:size(all_positions, 1)
    y_pos = all_positions(i, 1);
    item_type = all_positions(i, 2);
    area = all_positions(i, 3);
    
    % Distance from bottom
    dist_from_bottom = y_pos - y_bottom;
    
    % Determine color and label based on type
    if item_type == 1
        % Rebar (non-prestressed steel) - green
        lineColor = [0, 0.5, 0];
        label = sprintf('d_s = %.2f in', dist_from_bottom);
    else
        % Tendon (prestressed steel) - magenta
        tendon_num = all_positions(i, 4);
        lineColor = [0.6, 0, 0.6];
        label = sprintf('d_p%d = %.2f in', tendon_num, dist_from_bottom);
    end
    
    % Offset x position for multiple items to avoid overlap
    x_pos = x_dim_rebar - (i-1) * offset * 0.5;
    
    % Draw dimension line from bottom to rebar/tendon location
    drawDimensionLineColored(y_bottom, y_pos, x_pos, 'v', label, ext, fontsize, lineColor);
    
    % Draw a short horizontal line pointing to the rebar/tendon
    plot([x_pos, dims.x_min - offset*0.3], [y_pos, y_pos], ':', 'Color', lineColor, 'LineWidth', 0.5);
end

end

function drawIBeamDimensions(v, section, dims, offset, ext, gap, fontsize)
% Draw detailed dimensions for I-beam section

%% HORIZONTAL DIMENSIONS (Widths)
% Level 1: Total width at bottom (below section)
y_dim = dims.y_min - offset;
drawDimensionLine(dims.x_min, dims.x_max, y_dim, 'h', ...
    sprintf('%.1f in', dims.total_width), ext, fontsize);

% Level 2: Bottom flange width (if different from total)
if abs(dims.bot_flange_width - dims.total_width) > 0.1
    y_dim2 = dims.y_min - offset * 1.8;
    x_bf_left = dims.bot_flange_x(1);
    x_bf_right = dims.bot_flange_x(2);
    drawDimensionLine(x_bf_left, x_bf_right, y_dim2, 'h', ...
        sprintf('b_f = %.1f in', x_bf_right - x_bf_left), ext, fontsize);
end

% Top flange width (above section)
y_dim_top = dims.y_max + offset;
drawDimensionLine(dims.top_flange_x(1), dims.top_flange_x(2), y_dim_top, 'h', ...
    sprintf('b_f = %.1f in', dims.top_flange_width), ext, fontsize);

% Web width - show inside section
web_x_center = mean(dims.web_x);
web_y_center = mean(dims.web_y);
y_web_dim = web_y_center;
drawDimensionLine(dims.web_x(1), dims.web_x(2), y_web_dim, 'h', ...
    sprintf('b_w = %.1f in', dims.web_width), ext*0.5, fontsize-1, 'inside');

%% VERTICAL DIMENSIONS (Heights)
% Total height (right side, outermost)
x_dim = dims.x_max + offset * 2;
drawDimensionLine(dims.y_min, dims.y_max, x_dim, 'v', ...
    sprintf('h = %.1f in', dims.total_height), ext, fontsize);

% Individual component heights (right side, inner)
x_dim_inner = dims.x_max + offset * 0.8;

% Bottom flange height
drawDimensionLine(dims.bot_flange_y(1), dims.bot_flange_y(2), x_dim_inner, 'v', ...
    sprintf('%.1f in', dims.bot_flange_height), ext, fontsize-1);

% Web height
drawDimensionLine(dims.web_y(1), dims.web_y(2), x_dim_inner, 'v', ...
    sprintf('%.1f in', dims.web_height), ext, fontsize-1);

% Top flange height  
drawDimensionLine(dims.top_flange_y(1), dims.top_flange_y(2), x_dim_inner, 'v', ...
    sprintf('%.1f in', dims.top_flange_height), ext, fontsize-1);

%% CENTROID DISTANCES (Left side)
x_dim_left = dims.x_min - offset * 0.8;

% y_b (bottom to centroid)
drawDimensionLine(dims.y_min, section.yc, x_dim_left, 'v', ...
    sprintf('y_b = %.2f in', section.yb), ext, fontsize, 'yb');

% y_t (centroid to top)
drawDimensionLine(section.yc, dims.y_max, x_dim_left, 'v', ...
    sprintf('y_t = %.2f in', section.yt), ext, fontsize, 'yt');

%% FLANGE OVERHANG DIMENSIONS
if dims.top_flange_width > dims.web_width
    overhang_left = dims.web_x(1) - dims.top_flange_x(1);
    overhang_right = dims.top_flange_x(2) - dims.web_x(2);
    
    y_oh = dims.top_flange_y(1) - offset * 0.3;
    
    if abs(overhang_left) > 0.1
        drawDimensionLine(dims.top_flange_x(1), dims.web_x(1), y_oh, 'h', ...
            sprintf('%.1f in', overhang_left), ext*0.5, fontsize-2);
    end
    
    if abs(overhang_right) > 0.1
        drawDimensionLine(dims.web_x(2), dims.top_flange_x(2), y_oh, 'h', ...
            sprintf('%.1f in', overhang_right), ext*0.5, fontsize-2);
    end
end

end

function drawRectangularDimensions(v, section, dims, offset, ext, gap, fontsize)
% Draw dimensions for rectangular section

% Width (bottom)
y_dim = dims.y_min - offset;
drawDimensionLine(dims.x_min, dims.x_max, y_dim, 'h', ...
    sprintf('b = %.1f in', dims.width), ext, fontsize);

% Height (right)
x_dim = dims.x_max + offset;
drawDimensionLine(dims.y_min, dims.y_max, x_dim, 'v', ...
    sprintf('h = %.1f in', dims.height), ext, fontsize);

% Centroid distances (left)
x_dim_left = dims.x_min - offset;
drawDimensionLine(dims.y_min, section.yc, x_dim_left, 'v', ...
    sprintf('y_b = %.2f in', section.yb), ext, fontsize, 'yb');
drawDimensionLine(section.yc, dims.y_max, x_dim_left, 'v', ...
    sprintf('y_t = %.2f in', section.yt), ext, fontsize, 'yt');

end

function drawGeneralDimensions(v, section, dims, offset, ext, gap, fontsize)
% Draw dimensions for general polygon section

% Overall width (bottom)
y_dim = dims.y_min - offset;
drawDimensionLine(dims.x_min, dims.x_max, y_dim, 'h', ...
    sprintf('%.1f in', dims.total_width), ext, fontsize);

% Overall height (right)
x_dim = dims.x_max + offset;
drawDimensionLine(dims.y_min, dims.y_max, x_dim, 'v', ...
    sprintf('%.1f in', dims.total_height), ext, fontsize);

% Centroid distances (left)
x_dim_left = dims.x_min - offset;
drawDimensionLine(dims.y_min, section.yc, x_dim_left, 'v', ...
    sprintf('y_b = %.2f in', section.yb), ext, fontsize, 'yb');
drawDimensionLine(section.yc, dims.y_max, x_dim_left, 'v', ...
    sprintf('y_t = %.2f in', section.yt), ext, fontsize, 'yt');

% Draw dimensions at each unique y-level
for i = 2:length(dims.y_levels)-1
    y_level = dims.y_levels(i);
    tol = 0.1;
    idx_near = abs(v(:,2) - y_level) < dims.total_height * 0.05;
    
    if any(idx_near)
        x_at_level = v(idx_near, 1);
        if length(unique(x_at_level)) >= 2
            x_left = min(x_at_level);
            x_right = max(x_at_level);
            width_at_level = x_right - x_left;
            
            if abs(width_at_level - dims.total_width) > dims.total_width * 0.1
                plot([x_left, x_right], [y_level, y_level], 'b:', 'LineWidth', 0.5);
                text((x_left + x_right)/2, y_level + 1, ...
                    sprintf('%.1f in', width_at_level), ...
                    'FontSize', fontsize-1, 'HorizontalAlignment', 'center', 'Color', 'b');
            end
        end
    end
end

end

function drawDimensionLine(p1, p2, pos, orientation, label, ext, fontsize, style)
% Draw a single dimension line with extension lines and label

if nargin < 8
    style = 'normal';
end

lineColor = 'k';
lineWidth = 0.8;

% Colors for centroid distances
if strcmp(style, 'yb')
    lineColor = [0, 0, 0.7];  % Dark blue
elseif strcmp(style, 'yt')
    lineColor = [0.7, 0, 0];  % Dark red
end

if strcmp(orientation, 'h')
    % Horizontal dimension (measuring width)
    x1 = p1; x2 = p2;
    y = pos;
    
    if strcmp(style, 'inside')
        % Draw inside the section
        plot([x1, x2], [y, y], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        tick = ext * 0.5;
        plot([x1, x1], [y-tick, y+tick], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        plot([x2, x2], [y-tick, y+tick], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        text((x1+x2)/2, y + tick*2, label, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', fontsize, 'Color', lineColor);
    else
        % Extension lines (vertical)
        plot([x1, x1], [y, y + ext], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        plot([x2, x2], [y, y + ext], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        
        % Dimension line with arrows
        plot([x1, x2], [y, y], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        
        % Arrow heads
        arrow_len = abs(x2 - x1) * 0.03;
        arrow_len = max(arrow_len, 1);
        
        plot([x1, x1 + arrow_len], [y, y + arrow_len*0.4], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        plot([x1, x1 + arrow_len], [y, y - arrow_len*0.4], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        plot([x2, x2 - arrow_len], [y, y + arrow_len*0.4], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        plot([x2, x2 - arrow_len], [y, y - arrow_len*0.4], '-', 'Color', lineColor, 'LineWidth', lineWidth);
        
        text((x1+x2)/2, y - ext*0.3, label, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', fontsize, 'Color', lineColor);
    end
    
else
    % Vertical dimension (measuring height)
    y1 = p1; y2 = p2;
    x = pos;
    
    % Extension lines (horizontal)
    plot([x, x - ext], [y1, y1], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x, x - ext], [y2, y2], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    
    % Dimension line
    plot([x, x], [y1, y2], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    
    % Arrow heads
    arrow_len = abs(y2 - y1) * 0.03;
    arrow_len = max(arrow_len, 0.5);
    
    plot([x - arrow_len*0.4, x], [y1 + arrow_len, y1], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x + arrow_len*0.4, x], [y1 + arrow_len, y1], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x - arrow_len*0.4, x], [y2 - arrow_len, y2], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x + arrow_len*0.4, x], [y2 - arrow_len, y2], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    
    text(x + ext*0.3, (y1+y2)/2, label, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', fontsize, 'Color', lineColor, 'Rotation', 0);
end

end

function drawDimensionLineColored(p1, p2, pos, orientation, label, ext, fontsize, lineColor)
% Draw a dimension line with specified color

lineWidth = 0.8;

if strcmp(orientation, 'v')
    % Vertical dimension
    y1 = p1; y2 = p2;
    x = pos;
    
    % Extension lines (horizontal)
    plot([x, x - ext], [y1, y1], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x, x - ext], [y2, y2], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    
    % Dimension line
    plot([x, x], [y1, y2], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    
    % Arrow heads
    arrow_len = abs(y2 - y1) * 0.03;
    arrow_len = max(arrow_len, 0.5);
    
    plot([x - arrow_len*0.4, x], [y1 + arrow_len, y1], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x + arrow_len*0.4, x], [y1 + arrow_len, y1], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x - arrow_len*0.4, x], [y2 - arrow_len, y2], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x + arrow_len*0.4, x], [y2 - arrow_len, y2], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    
    text(x + ext*0.3, (y1+y2)/2, label, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', fontsize, 'Color', lineColor, 'Rotation', 0);
else
    % Horizontal dimension
    x1 = p1; x2 = p2;
    y = pos;
    
    plot([x1, x1], [y, y + ext], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x2, x2], [y, y + ext], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x1, x2], [y, y], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    
    arrow_len = abs(x2 - x1) * 0.03;
    arrow_len = max(arrow_len, 1);
    
    plot([x1, x1 + arrow_len], [y, y + arrow_len*0.4], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x1, x1 + arrow_len], [y, y - arrow_len*0.4], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x2, x2 - arrow_len], [y, y + arrow_len*0.4], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot([x2, x2 - arrow_len], [y, y - arrow_len*0.4], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    
    text((x1+x2)/2, y - ext*0.3, label, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', fontsize, 'Color', lineColor);
end

end
