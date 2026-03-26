function [beam, section, materials, prestress, reinforcement, loads] = quickInputPE(P, e, L, section_vertices, varargin)
% QUICKINPUTPE - Quick input for prestressed beam analysis using P and e
%
% Simplified input function where prestress is defined by:
%   P - Total effective prestress force (kips)
%   e - Eccentricity profile or value (in, positive below centroid)
%
% Syntax:
%   [beam, section, materials, prestress, reinforcement, loads] = quickInputPE(P, e, L, section_vertices)
%   [...] = quickInputPE(P, e, L, section_vertices, 'Name', Value, ...)
%
% Required Inputs:
%   P               - Effective prestress force (kips) - scalar or [1 x n] array
%   e               - Eccentricity (in) - scalar, [3x1] for [e_start; e_mid; e_end], or [1 x n] array
%   L               - Beam length (in)
%   section_vertices - [n x 2] matrix of section coordinates (counterclockwise)
%
% Optional Name-Value Pairs:
%   'ProfileType'       - 'linear', 'parabolic', 'harped', 'custom' (default: 'parabolic')
%   'DrapePoint'        - Drape point as fraction of span for harped profiles (default: 0.4)
%   'BondingType'       - 'full', 'unbonded', 'partial' (default: 'full')
%   'BondedZones'       - [n x 2] matrix of [start, end] fractions for partial bonding
%   'Losses'            - Prestress losses as fraction (default: 0.15)
%   'fc'                - Concrete strength in ksi (default: 6)
%   'DistributedLoad'   - Uniform distributed load in kip/in (default: 0)
%   'PointLoads'        - [n x 3] matrix of [x, P, M] point loads
%   'SupportType'       - 'simple', 'cantilever' (default: 'simple')
%   'SupportLocations'  - [x_left, x_right] for simple supports (default: [0, L])
%   'NumSegments'       - Number of analysis segments (default: 100)
%   'RebarTop'          - [n x 3] matrix of [x, y, area] for top steel
%   'RebarBottom'       - [n x 3] matrix of [x, y, area] for bottom steel
%
% Examples:
%   % Simple rectangular beam with constant eccentricity
%   rect = [0 0; 12 0; 12 24; 0 24];
%   [beam, section, materials, prestress, reinforcement, loads] = ...
%       quickInputPE(200, 6, 360, rect);
%
%   % I-section with parabolic tendon profile
%   I_section = [0 0; 24 0; 24 6; 15 6; 15 30; 24 30; 24 36; 0 36; 0 30; 9 30; 9 6; 0 6];
%   [beam, section, materials, prestress, reinforcement, loads] = ...
%       quickInputPE(350, [4; 12; 4], 480, I_section, 'ProfileType', 'parabolic');
%
%   % With partial bonding
%   [beam, section, materials, prestress, reinforcement, loads] = ...
%       quickInputPE(350, [4; 12; 4], 480, I_section, ...
%       'BondingType', 'partial', 'BondedZones', [0, 0.2; 0.4, 0.6; 0.8, 1.0]);

%% Parse inputs
p = inputParser;
addRequired(p, 'P', @isnumeric);
addRequired(p, 'e', @isnumeric);
addRequired(p, 'L', @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'section_vertices', @(x) isnumeric(x) && size(x,2) == 2);

addParameter(p, 'ProfileType', 'parabolic', @ischar);
addParameter(p, 'DrapePoint', 0.4, @isnumeric);
addParameter(p, 'BondingType', 'full', @ischar);
addParameter(p, 'BondedZones', [0, 1], @isnumeric);
addParameter(p, 'Losses', 0.15, @isnumeric);
addParameter(p, 'fc', 6, @isnumeric);
addParameter(p, 'DistributedLoad', 0, @isnumeric);
addParameter(p, 'PointLoads', [], @isnumeric);
addParameter(p, 'SupportType', 'simple', @ischar);
addParameter(p, 'SupportLocations', [], @isnumeric);
addParameter(p, 'NumSegments', 100, @isnumeric);
addParameter(p, 'RebarTop', [], @isnumeric);
addParameter(p, 'RebarBottom', [], @isnumeric);

parse(p, P, e, L, section_vertices, varargin{:});

opts = p.Results;

%% Setup beam geometry
beam.L = L;
beam.num_segments = opts.NumSegments;
beam.x = linspace(0, L, beam.num_segments + 1);

%% Calculate section properties
section = createPolygonalSection(section_vertices, 'Name', 'User Section');

%% Setup materials
materials.fc = opts.fc;
materials.Ec = 57 * sqrt(materials.fc * 1000);
materials.fci = 0.75 * materials.fc;
materials.fr = 7.5 * sqrt(materials.fc * 1000) / 1000;

materials.fpu = 270;
materials.fpy = 0.9 * materials.fpu;
materials.Eps = 28500;

materials.fy = 60;
materials.Es = 29000;

%% Setup prestressing
% Convert P to effective force (assuming P is already effective after losses)
% If user wants to specify initial prestress, they should account for losses

% Determine eccentricity profile
if isscalar(e)
    % Constant eccentricity
    e_start = e;
    e_mid = e;
    e_end = e;
elseif length(e) == 3
    % [e_start; e_mid; e_end]
    e_start = e(1);
    e_mid = e(2);
    e_end = e(3);
else
    % Custom profile - must match x array length
    if length(e) ~= length(beam.x)
        error('Custom eccentricity array must have same length as x array (%d)', length(beam.x));
    end
end

% Create single equivalent tendon
tendon.Aps = P / (materials.fpu * (1 - opts.Losses) * 0.75);  % Back-calculate Aps
tendon.fpi = 0.75 * materials.fpu;  % 75% of fpu
tendon.profile_type = opts.ProfileType;

if ~isscalar(e) && length(e) == length(beam.x)
    tendon.profile_type = 'custom';
    tendon.x_profile = beam.x;
    tendon.e_profile = e;
    tendon.e_start = e(1);
    tendon.e_mid = e(round(end/2));
    tendon.e_end = e(end);
else
    tendon.e_start = e_start;
    tendon.e_mid = e_mid;
    tendon.e_end = e_end;
end

tendon.drape_point = opts.DrapePoint;

% Bonding
tendon.bonding.type = opts.BondingType;
tendon.bonding.bonded_zones = opts.BondedZones;

prestress.tendons = {tendon};
prestress.losses = opts.Losses;

% Process tendon profile
prestress = processTendonProfiles(prestress, beam, section);

% Override with direct P if provided as array
if ~isscalar(P)
    if length(P) == length(beam.x)
        % Use custom P distribution
        prestress.tendons{1}.P_effective = P;
    else
        error('Custom P array must have same length as x array (%d)', length(beam.x));
    end
end

%% Setup reinforcement
reinforcement.longitudinal = {};

if ~isempty(opts.RebarTop)
    for i = 1:size(opts.RebarTop, 1)
        rebar.x_positions = opts.RebarTop(i, 1);
        rebar.y_position = opts.RebarTop(i, 2);
        rebar.areas = opts.RebarTop(i, 3);
        rebar.type = 'compression';
        reinforcement.longitudinal{end+1} = rebar;
    end
end

if ~isempty(opts.RebarBottom)
    for i = 1:size(opts.RebarBottom, 1)
        rebar.x_positions = opts.RebarBottom(i, 1);
        rebar.y_position = opts.RebarBottom(i, 2);
        rebar.areas = opts.RebarBottom(i, 3);
        rebar.type = 'tension';
        reinforcement.longitudinal{end+1} = rebar;
    end
end

reinforcement.stirrups.Av = 0.22;
reinforcement.stirrups.spacing = 6;
reinforcement.stirrups.fy = 60;

%% Setup loads
loads.self_weight = [];
loads.concrete_density = 150/1728;

% Distributed loads
if opts.DistributedLoad ~= 0
    loads.distributed = [0, L, opts.DistributedLoad, opts.DistributedLoad];
else
    loads.distributed = [];
end

% Point loads
if ~isempty(opts.PointLoads)
    loads.point = opts.PointLoads;
else
    loads.point = [];
end

% Supports
loads.support_type = opts.SupportType;
if isempty(opts.SupportLocations)
    loads.supports = [0, L];
else
    loads.supports = opts.SupportLocations;
end

%% Display summary
fprintf('\n========================================\n');
fprintf('   QUICK INPUT (P, e) SUMMARY\n');
fprintf('========================================\n');
fprintf('Prestress Force P: %.1f kips\n', P(1));
if isscalar(e)
    fprintf('Eccentricity e: %.2f in (constant)\n', e);
else
    fprintf('Eccentricity e: %.2f to %.2f in (%s profile)\n', ...
        min(e), max(e), opts.ProfileType);
end
fprintf('Beam Length L: %.1f in (%.2f ft)\n', L, L/12);
fprintf('Section Area: %.1f in²\n', section.A);
fprintf('Bonding: %s\n', opts.BondingType);
fprintf('========================================\n\n');

end

%% HELPER FUNCTION
function prestress = processTendonProfiles(prestress, beam, section)
x = beam.x;
L = beam.L;
yc = section.yc;

for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};
    
    e_s = tendon.e_start;
    e_m = tendon.e_mid;
    e_e = tendon.e_end;
    
    switch tendon.profile_type
        case 'linear'
            e = e_s + (e_e - e_s) * x / L;
        case 'parabolic'
            a = 2 * (e_s - 2*e_m + e_e) / L^2;
            b = (e_e - e_s) / L - a * L;
            c = e_s;
            e = a * x.^2 + b * x + c;
        case 'harped'
            x_drape = tendon.drape_point * L;
            e = zeros(size(x));
            mask1 = x <= x_drape;
            mask2 = x > x_drape;
            e(mask1) = e_s + (e_m - e_s) * x(mask1) / x_drape;
            e(mask2) = e_m + (e_e - e_m) * (x(mask2) - x_drape) / (L - x_drape);
        case 'custom'
            e = interp1(tendon.x_profile, tendon.e_profile, x, 'linear');
        otherwise
            error('Unknown profile type');
    end
    
    prestress.tendons{i}.e = e;
    prestress.tendons{i}.y = yc - e;
    prestress.tendons{i}.bonded = generateBondingMask(tendon.bonding, x, L);
end
end

function bonded = generateBondingMask(bonding, x, L)
bonded = false(size(x));
switch bonding.type
    case 'full'
        bonded(:) = true;
    case 'unbonded'
        bonded(:) = false;
    case 'partial'
        for i = 1:size(bonding.bonded_zones, 1)
            x_start = bonding.bonded_zones(i, 1) * L;
            x_end = bonding.bonded_zones(i, 2) * L;
            bonded = bonded | (x >= x_start & x <= x_end);
        end
end
end
