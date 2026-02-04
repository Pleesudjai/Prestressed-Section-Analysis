function [results] = analyzePrestressedBeam(beam, section, materials, prestress, reinforcement, loads)
% ANALYZEPRESTRESSEDBEAM - Analyze prestressed concrete beam
% Calculates axial force (N), shear force (V), and bending moment (M) diagrams
%
% Inputs:
%   beam - Beam geometry structure
%   section - Cross-section properties structure
%   materials - Material properties structure
%   prestress - Prestressing tendon information structure
%   reinforcement - Non-prestressed reinforcement structure
%   loads - Applied loads structure
%
% Output:
%   results - Structure containing all analysis results

%% Initialize results structure
x = beam.x;
n = length(x);
L = beam.L;

results.x = x;
results.L = L;

% Initialize force arrays
results.N = zeros(1, n);        % Axial force (kips)
results.V = zeros(1, n);        % Shear force (kips)
results.M = zeros(1, n);        % Bending moment (kip-in)

% Prestress effects
results.P = zeros(1, n);        % Effective prestress force (kips)
results.e = zeros(1, n);        % Effective eccentricity (in)
results.M_prestress = zeros(1, n);  % Moment due to prestress (kip-in)

% Component breakdown
results.V_dead = zeros(1, n);
results.V_live = zeros(1, n);
results.V_prestress = zeros(1, n);
results.M_dead = zeros(1, n);
results.M_live = zeros(1, n);

%% Calculate self-weight if not provided
if isempty(loads.self_weight)
    loads.self_weight = section.A * loads.concrete_density;
end

%% Calculate prestress effects
[results.P, results.e, results.M_prestress, results.V_prestress] = ...
    calculatePrestressEffects(beam, section, prestress);

%% Calculate reactions and internal forces
[results.V, results.M, results.V_dead, results.V_live, results.M_dead, results.M_live, results.reactions] = ...
    calculateInternalForces(beam, loads, section);

%% Combine effects - Axial force from prestress
% For simply supported beam, external axial force is typically zero
% Prestress creates internal compression
results.N = -results.P;  % Compression negative in structural convention

%% Calculate stresses at critical sections
results.stresses = calculateStresses(results, section, materials);

%% Calculate capacity envelopes (optional detailed analysis)
results.capacity = calculateCapacity(section, materials, prestress, reinforcement);

%% Store additional information
results.section = section;
results.materials = materials;
results.prestress = prestress;
results.loads = loads;

end

%% PRESTRESS EFFECTS
function [P, e_eff, M_prestress, V_prestress] = calculatePrestressEffects(beam, section, prestress)
% Calculate effective prestress force, eccentricity, and equivalent loads

x = beam.x;
n = length(x);
L = beam.L;

P = zeros(1, n);
e_eff = zeros(1, n);
M_prestress = zeros(1, n);
V_prestress = zeros(1, n);

% Sum contributions from all tendons
for i = 1:length(prestress.tendons)
    tendon = prestress.tendons{i};
    
    % Effective prestress force (after losses)
    Pe = tendon.Aps * tendon.fpi * (1 - prestress.losses);
    
    % Eccentricity profile
    e = tendon.e;
    
    % Add to totals (weighted by bonding for moment effects)
    P = P + Pe * ones(1, n);  % Axial force always acts
    e_eff = e_eff + Pe * e;   % Weighted eccentricity
    
    % Primary moment from prestress: M = P * e
    M_primary = Pe * e;
    
    % Calculate equivalent loads for secondary effects (if applicable)
    % For simple spans, secondary moments are zero
    % For continuous beams, would need additional analysis
    
    M_prestress = M_prestress + M_primary;
    
    % Shear from prestress (due to tendon inclination)
    % V = P * de/dx
    de_dx = gradient(e, x);
    V_prestress = V_prestress + Pe * de_dx;
end

% Normalize effective eccentricity by total force
total_P = sum(cellfun(@(t) t.Aps * t.fpi * (1 - prestress.losses), prestress.tendons));
if total_P > 0
    e_eff = e_eff / total_P;
end

end

%% INTERNAL FORCES
function [V, M, V_dead, V_live, M_dead, M_live, reactions] = calculateInternalForces(beam, loads, section)
% Calculate shear and moment diagrams from applied loads

x = beam.x;
n = length(x);
L = beam.L;
dx = x(2) - x(1);

% Initialize arrays
V = zeros(1, n);
M = zeros(1, n);
V_dead = zeros(1, n);
V_live = zeros(1, n);
M_dead = zeros(1, n);
M_live = zeros(1, n);

%% Process distributed loads
w = zeros(1, n);  % Total distributed load at each point
w_dead = zeros(1, n);
w_live = zeros(1, n);

% Add self-weight (dead load)
w_dead = w_dead - loads.self_weight * ones(1, n);

% Add user-defined distributed loads
for i = 1:size(loads.distributed, 1)
    x_start = loads.distributed(i, 1);
    x_end = loads.distributed(i, 2);
    w_start = loads.distributed(i, 3);
    w_end = loads.distributed(i, 4);
    
    % Interpolate load intensity along the segment
    mask = (x >= x_start) & (x <= x_end);
    if any(mask)
        w_segment = interp1([x_start, x_end], [w_start, w_end], x(mask), 'linear');
        
        % Classify as dead or live based on sign or position
        % For simplicity, all distributed loads go to dead load category
        w_dead(mask) = w_dead(mask) + w_segment;
    end
end

w = w_dead + w_live;

%% Calculate reactions based on support type
switch loads.support_type
    case 'simple'
        % Simple supports at specified locations
        x_left = loads.supports(1);
        x_right = loads.supports(2);
        span = x_right - x_left;
        
        % Calculate total load and moment about left support
        total_load = trapz(x, w);
        moment_left = trapz(x, w .* (x - x_left));
        
        % Add point loads
        for i = 1:size(loads.point, 1)
            x_p = loads.point(i, 1);
            P_p = loads.point(i, 2);
            M_p = loads.point(i, 3);
            
            total_load = total_load + P_p;
            moment_left = moment_left + P_p * (x_p - x_left) + M_p;
        end
        
        % Reactions
        R_right = moment_left / span;
        R_left = -total_load - R_right;
        
        reactions.left = R_left;
        reactions.right = R_right;
        reactions.x_left = x_left;
        reactions.x_right = x_right;
        
    case 'cantilever'
        % Cantilever from left support
        x_support = loads.supports(1);
        
        % Fixed-end reactions
        total_load = trapz(x, w);
        moment_support = trapz(x, w .* (x - x_support));
        
        for i = 1:size(loads.point, 1)
            x_p = loads.point(i, 1);
            P_p = loads.point(i, 2);
            M_p = loads.point(i, 3);
            
            total_load = total_load + P_p;
            moment_support = moment_support + P_p * (x_p - x_support) + M_p;
        end
        
        reactions.vertical = -total_load;
        reactions.moment = -moment_support;
        reactions.x_support = x_support;
        
    otherwise
        error('Unsupported support type: %s', loads.support_type);
end

%% Calculate shear and moment using integration
% Using the relationship: dV/dx = -w, dM/dx = V

switch loads.support_type
    case 'simple'
        % Start from left support
        V(1) = reactions.left;
        M(1) = 0;
        
        % Find indices of support locations
        [~, idx_left] = min(abs(x - reactions.x_left));
        [~, idx_right] = min(abs(x - reactions.x_right));
        
        % Integrate shear to get moment
        for i = 2:n
            % Check for point loads at this location
            P_point = 0;
            M_point = 0;
            for j = 1:size(loads.point, 1)
                if abs(x(i) - loads.point(j, 1)) < dx/2
                    P_point = P_point + loads.point(j, 2);
                    M_point = M_point + loads.point(j, 3);
                end
            end
            
            % Check for support reactions
            if abs(x(i) - reactions.x_left) < dx/2
                P_point = P_point + reactions.left;
            end
            if abs(x(i) - reactions.x_right) < dx/2
                P_point = P_point + reactions.right;
            end
            
            % Integrate: V(i) = V(i-1) - integral(w)
            V(i) = V(i-1) - trapz(x(i-1:i), w(i-1:i)) + P_point;
            
            % Integrate: M(i) = M(i-1) + integral(V)
            M(i) = M(i-1) + trapz(x(i-1:i), V(i-1:i)) + M_point;
        end
        
    case 'cantilever'
        % Start from support
        [~, idx_support] = min(abs(x - reactions.x_support));
        
        V(idx_support) = reactions.vertical;
        M(idx_support) = reactions.moment;
        
        % Integrate forward
        for i = idx_support+1:n
            P_point = 0;
            M_point = 0;
            for j = 1:size(loads.point, 1)
                if abs(x(i) - loads.point(j, 1)) < dx/2
                    P_point = P_point + loads.point(j, 2);
                    M_point = M_point + loads.point(j, 3);
                end
            end
            
            V(i) = V(i-1) - trapz(x(i-1:i), w(i-1:i)) + P_point;
            M(i) = M(i-1) + trapz(x(i-1:i), V(i-1:i)) + M_point;
        end
        
        % Integrate backward from support
        for i = idx_support-1:-1:1
            P_point = 0;
            M_point = 0;
            for j = 1:size(loads.point, 1)
                if abs(x(i+1) - loads.point(j, 1)) < dx/2
                    P_point = P_point + loads.point(j, 2);
                    M_point = M_point + loads.point(j, 3);
                end
            end
            
            V(i) = V(i+1) + trapz(x(i:i+1), w(i:i+1)) - P_point;
            M(i) = M(i+1) - trapz(x(i:i+1), V(i:i+1)) - M_point;
        end
end

% Separate dead and live load effects (simplified approach)
% In practice, would need separate calculations
V_dead = V;  % Placeholder - all loads treated as dead for now
M_dead = M;
V_live = zeros(1, n);
M_live = zeros(1, n);

end

%% STRESS CALCULATIONS
function stresses = calculateStresses(results, section, materials)
% Calculate stresses at critical locations

n = length(results.x);

% Initialize stress arrays
stresses.f_top = zeros(1, n);      % Stress at top fiber (ksi)
stresses.f_bot = zeros(1, n);      % Stress at bottom fiber (ksi)
stresses.f_top_prestress = zeros(1, n);
stresses.f_bot_prestress = zeros(1, n);
stresses.f_top_total = zeros(1, n);
stresses.f_bot_total = zeros(1, n);

A = section.A;
Ix = section.Ix;
yt = section.yt;
yb = section.yb;

for i = 1:n
    P = results.P(i);           % Prestress force (kips)
    e = results.e(i);           % Eccentricity (in)
    M = results.M(i);           % External moment (kip-in)
    
    % Prestress stresses: f = -P/A ± P*e*y/I
    % Convention: compression negative
    f_prestress_top = -P/A + P*e*yt/Ix;    % Top fiber
    f_prestress_bot = -P/A - P*e*yb/Ix;    % Bottom fiber
    
    % External load stresses: f = M*y/I
    % For positive moment (sagging): tension at bottom
    f_load_top = -M*yt/Ix;   % Compression at top for positive M
    f_load_bot = M*yb/Ix;    % Tension at bottom for positive M
    
    % Store results
    stresses.f_top_prestress(i) = f_prestress_top;
    stresses.f_bot_prestress(i) = f_prestress_bot;
    stresses.f_top(i) = f_load_top;
    stresses.f_bot(i) = f_load_bot;
    stresses.f_top_total(i) = f_prestress_top + f_load_top;
    stresses.f_bot_total(i) = f_prestress_bot + f_load_bot;
end

% Allowable stresses
stresses.fc_allow_compression = -0.45 * materials.fc;  % Compression limit
stresses.fc_allow_tension = materials.fr;              % Tension limit

end

%% CAPACITY CALCULATIONS
function capacity = calculateCapacity(section, materials, prestress, reinforcement)
% Calculate moment capacity at critical sections

% Nominal moment capacity (simplified)
% Using rectangular stress block method

fc = materials.fc;
fy = materials.fy;
fpu = materials.fpu;
fpy = materials.fpy;  % Yield strength of prestressing steel
beta1 = 0.85 - 0.05 * (fc - 4);
beta1 = max(0.65, min(0.85, beta1));

% Total prestressing steel area
Aps_total = sum(cellfun(@(t) t.Aps, prestress.tendons));

% Average tendon depth (approximate)
yc = section.yc;
e_avg = mean(cellfun(@(t) mean(t.e), prestress.tendons));
dp = yc + e_avg;

% Effective prestress
fpe = prestress.tendons{1}.fpi * (1 - prestress.losses);

% Stress in prestressing steel at nominal strength (approximate)
fps = fpu * (1 - 0.28 * (fpu / fc) * (Aps_total / (section.A * dp / section.yb)));
fps = min(fps, fpy);

% Calculate c (neutral axis depth) - iterative or approximate
% Simplified: assume fps ≈ fpe + 15 ksi for bonded tendons
fps_approx = min(fpe + 15, fpu);

% Get width at different depths (approximate as average width)
b = section.A / (section.yt + section.yb);

% Equilibrium: 0.85*fc*beta1*c*b = Aps*fps
c = Aps_total * fps_approx / (0.85 * fc * beta1 * b);
a = beta1 * c;

% Nominal moment capacity
Mn = Aps_total * fps_approx * (dp - a/2);

% Phi factor for prestressed members
epsilon_t = 0.003 * (dp - c) / c;  % Tensile strain
if epsilon_t >= 0.005
    phi = 0.9;  % Tension-controlled
elseif epsilon_t <= 0.002
    phi = 0.65; % Compression-controlled
else
    phi = 0.65 + (epsilon_t - 0.002) * (0.9 - 0.65) / (0.005 - 0.002);
end

capacity.Mn = Mn;
capacity.phi = phi;
capacity.phi_Mn = phi * Mn;
capacity.c = c;
capacity.a = a;
capacity.fps = fps_approx;
capacity.beta1 = beta1;

end
