function section = createPolygonalSection(vertices, varargin)
% CREATEPOLYGONALSECTION - Create section properties from polygon vertices
%
% Universal function to define any polygonal cross-section
%
% Syntax:
%   section = createPolygonalSection(vertices)
%   section = createPolygonalSection(vertices, 'PropertyName', PropertyValue, ...)
%
% Inputs:
%   vertices - [n x 2] matrix of [x, y] coordinates (counterclockwise)
%
% Optional Name-Value Pairs:
%   'Holes' - Cell array of hole vertex matrices (clockwise orientation)
%   'Name'  - String name for the section
%
% Output:
%   section - Structure with calculated section properties
%
% Examples:
%   % Rectangle 12" x 24"
%   section = createPolygonalSection([0 0; 12 0; 12 24; 0 24]);
%
%   % I-Section
%   vertices = [0 0; 24 0; 24 6; 15 6; 15 30; 24 30; 24 36; 0 36; 0 30; 9 30; 9 6; 0 6];
%   section = createPolygonalSection(vertices);
%
%   % Box section (with hole)
%   outer = [0 0; 24 0; 24 24; 0 24];
%   inner = [4 4; 4 20; 20 20; 20 4];  % Clockwise for hole
%   section = createPolygonalSection(outer, 'Holes', {inner});

%% Parse inputs
p = inputParser;
addRequired(p, 'vertices', @(x) isnumeric(x) && size(x,2) == 2);
addParameter(p, 'Holes', {}, @iscell);
addParameter(p, 'Name', 'Custom Section', @ischar);
parse(p, vertices, varargin{:});

holes = p.Results.Holes;
section.name = p.Results.Name;

%% Store original vertices
section.vertices = vertices;
section.holes = holes;

%% Ensure counterclockwise orientation for outer boundary
if ~isCounterclockwise(vertices)
    vertices = flipud(vertices);
    section.vertices = vertices;
    warning('Vertices were reordered to counterclockwise orientation');
end

%% Calculate properties for outer boundary
[A_outer, Cx_outer, Cy_outer, Ix_outer, Iy_outer, ~] = ...
    calculatePolygonProperties(vertices);

%% Calculate properties for holes (subtract)
A_holes = 0;
Sx_holes = 0;
Sy_holes = 0;

for i = 1:length(holes)
    hole = holes{i};
    
    % Ensure clockwise orientation for holes
    if isCounterclockwise(hole)
        hole = flipud(hole);
        section.holes{i} = hole;
    end
    
    [A_h, Cx_h, Cy_h, ~, ~, ~] = calculatePolygonProperties(hole);
    A_h = abs(A_h);
    
    A_holes = A_holes + A_h;
    Sx_holes = Sx_holes + A_h * Cx_h;
    Sy_holes = Sy_holes + A_h * Cy_h;
end

%% Net section properties
section.A = A_outer - A_holes;

if section.A <= 0
    error('Net section area is non-positive. Check hole definitions.');
end

% Net centroid
Sx_outer = A_outer * Cx_outer;
Sy_outer = A_outer * Cy_outer;

section.xc = (Sx_outer - Sx_holes) / section.A;
section.yc = (Sy_outer - Sy_holes) / section.A;

%% Calculate net moment of inertia about the net centroid
[section.Ix, section.Iy, section.Ixy] = calculateNetInertia(vertices, holes, section.xc, section.yc);

%% Calculate extreme fiber distances
y_all = vertices(:, 2);
x_all = vertices(:, 1);

y_max = max(y_all);
y_min = min(y_all);
x_max = max(x_all);
x_min = min(x_all);

section.yt = y_max - section.yc;
section.yb = section.yc - y_min;
section.xr = x_max - section.xc;
section.xl = section.xc - x_min;

%% Section moduli
section.St = section.Ix / section.yt;
section.Sb = section.Ix / section.yb;
section.Sl = section.Iy / section.xl;
section.Sr = section.Iy / section.xr;

%% Radius of gyration
section.rx = sqrt(section.Ix / section.A);
section.ry = sqrt(section.Iy / section.A);

%% Kern dimensions (for prestress)
section.kt = section.Sb / section.A;
section.kb = section.St / section.A;

%% Display summary
fprintf('\n--- Section Properties: %s ---\n', section.name);
fprintf('Area (A):              %.2f in²\n', section.A);
fprintf('Centroid (xc, yc):     (%.2f, %.2f) in\n', section.xc, section.yc);
fprintf('Moment of Inertia Ix:  %.0f in⁴\n', section.Ix);
fprintf('Moment of Inertia Iy:  %.0f in⁴\n', section.Iy);
fprintf('Distance to top (yt):  %.2f in\n', section.yt);
fprintf('Distance to bot (yb):  %.2f in\n', section.yb);
fprintf('Section Modulus St:    %.0f in³\n', section.St);
fprintf('Section Modulus Sb:    %.0f in³\n', section.Sb);
fprintf('Kern kt, kb:           %.2f, %.2f in\n', section.kt, section.kb);
fprintf('-----------------------------------\n\n');

end

%% HELPER FUNCTIONS

function ccw = isCounterclockwise(vertices)
n = size(vertices, 1);
v = [vertices; vertices(1,:)];
signed_area = 0;
for i = 1:n
    signed_area = signed_area + (v(i,1) * v(i+1,2) - v(i+1,1) * v(i,2));
end
ccw = signed_area > 0;
end

function [A, Cx, Cy, Ix, Iy, Ixy] = calculatePolygonProperties(vertices)
n = size(vertices, 1);
v = [vertices; vertices(1,:)];

A = 0;
for i = 1:n
    A = A + (v(i,1) * v(i+1,2) - v(i+1,1) * v(i,2));
end
A = A / 2;

Cx = 0; Cy = 0;
for i = 1:n
    factor = (v(i,1) * v(i+1,2) - v(i+1,1) * v(i,2));
    Cx = Cx + (v(i,1) + v(i+1,1)) * factor;
    Cy = Cy + (v(i,2) + v(i+1,2)) * factor;
end
Cx = Cx / (6 * A);
Cy = Cy / (6 * A);

Ix = 0; Iy = 0; Ixy = 0;
for i = 1:n
    xi = v(i,1); yi = v(i,2);
    xi1 = v(i+1,1); yi1 = v(i+1,2);
    cross = (xi * yi1 - xi1 * yi);
    Ix = Ix + (yi^2 + yi*yi1 + yi1^2) * cross;
    Iy = Iy + (xi^2 + xi*xi1 + xi1^2) * cross;
    Ixy = Ixy + (xi*yi1 + 2*xi*yi + 2*xi1*yi1 + xi1*yi) * cross;
end

Ix = Ix / 12;
Iy = Iy / 12;
Ixy = Ixy / 24;

Ix = abs(Ix) - abs(A) * Cy^2;
Iy = abs(Iy) - abs(A) * Cx^2;
Ixy = Ixy - A * Cx * Cy;
end

function [Ix, Iy, Ixy] = calculateNetInertia(outer, holes, xc, yc)
[A_o, Cx_o, Cy_o, Ix_o, Iy_o, Ixy_o] = calculatePolygonPropertiesRaw(outer);

dx_o = Cx_o - xc;
dy_o = Cy_o - yc;
Ix = Ix_o + abs(A_o) * dy_o^2;
Iy = Iy_o + abs(A_o) * dx_o^2;
Ixy = Ixy_o + A_o * dx_o * dy_o;

for i = 1:length(holes)
    [A_h, Cx_h, Cy_h, Ix_h, Iy_h, Ixy_h] = calculatePolygonPropertiesRaw(holes{i});
    A_h = abs(A_h);
    Ix_h = abs(Ix_h);
    Iy_h = abs(Iy_h);
    
    dx_h = Cx_h - xc;
    dy_h = Cy_h - yc;
    
    Ix = Ix - (Ix_h + A_h * dy_h^2);
    Iy = Iy - (Iy_h + A_h * dx_h^2);
    Ixy = Ixy - (Ixy_h + A_h * dx_h * dy_h);
end
end

function [A, Cx, Cy, Ix, Iy, Ixy] = calculatePolygonPropertiesRaw(vertices)
n = size(vertices, 1);
v = [vertices; vertices(1,:)];

A = 0;
for i = 1:n
    A = A + (v(i,1) * v(i+1,2) - v(i+1,1) * v(i,2));
end
A = A / 2;

Cx = 0; Cy = 0;
for i = 1:n
    factor = (v(i,1) * v(i+1,2) - v(i+1,1) * v(i,2));
    Cx = Cx + (v(i,1) + v(i+1,1)) * factor;
    Cy = Cy + (v(i,2) + v(i+1,2)) * factor;
end
Cx = Cx / (6 * A);
Cy = Cy / (6 * A);

Ix = 0; Iy = 0; Ixy = 0;
for i = 1:n
    xi = v(i,1); yi = v(i,2);
    xi1 = v(i+1,1); yi1 = v(i+1,2);
    cross = (xi * yi1 - xi1 * yi);
    Ix = Ix + (yi^2 + yi*yi1 + yi1^2) * cross;
    Iy = Iy + (xi^2 + xi*xi1 + xi1^2) * cross;
    Ixy = Ixy + (xi*yi1 + 2*xi*yi + 2*xi1*yi1 + xi1*yi) * cross;
end

Ix = Ix / 12;
Iy = Iy / 12;
Ixy = Ixy / 24;

Ix = abs(Ix) - abs(A) * Cy^2;
Iy = abs(Iy) - abs(A) * Cx^2;
Ixy = Ixy - A * Cx * Cy;
end
