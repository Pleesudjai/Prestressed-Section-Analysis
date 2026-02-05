function [results] = analyzePrestressedBeam(beam, section, materials, prestress, reinforcement, loads)
% ANALYZEPRESTRESSEDBEAM - Analyze prestressed concrete beam
% Calculates axial force (N), shear force (V), and bending moment (M) diagrams
%
% SIGN CONVENTION (Universal Beam Convention):
%   Positive Moment (+M) = Sagging (tension at bottom, compression at top)
%   Negative Moment (-M) = Hogging (tension at top, compression at bottom)
%   Positive eccentricity (e > 0) = tendon BELOW centroid
%   Prestress force P > 0 (compression)
%
% Prestress Moment:
%   Tendon below centroid (e > 0) creates HOGGING → M_prestress = -P*e
%   Tendon above centroid (e < 0) creates SAGGING → M_prestress = -P*e (still works)
%
% Boundary Conditions for Simply Supported Beam:
%   At x = 0 (left support):
%     - M(0) = 0 (moment boundary condition)
%     - V(0+) = R_left (shear jumps by reaction)
%   At x = L (right support):
%     - M(L) = 0 (moment boundary condition, satisfied by correct reactions)
%     - V(L-) = -R_right (shear before support)
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

%% PRESTRESS EFFECTS - CORRECTED SIGN
function [P, e_eff, M_prestress, V_prestress] = calculatePrestressEffects(beam, section, prestress)
% Calculate effective prestress force, eccentricity, and equivalent loads
%
% SIGN CONVENTION:
%   e > 0 means tendon is BELOW centroid
%   P > 0 (always positive, it's a force magnitude)
%   M_prestress = -P * e  (negative = hogging when tendon below centroid)
%
% Physical reasoning:
%   Tendon below centroid pulls the bottom fiber into compression
%   more than the top → beam cambers upward → HOGGING moment → NEGATIVE M

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
   
    % Eccentricity profile (positive = below centroid)
    e = tendon.e;
    
    % Add to totals
    P = P + Pe * ones(1, n);  % Axial force always acts
    e_eff = e_eff + Pe * e;   % Weighted eccentricity
    
    % -------------------------------------------------------
    % CORRECTED: Primary moment from prestress
    %   M_prestress = -P * e
    %   When e > 0 (tendon below CG): M < 0 → hogging (camber)
    %   When e < 0 (tendon above CG): M > 0 → sagging
    % -------------------------------------------------------
    M_primary = -Pe * e;
    
    M_prestress = M_prestress + M_primary;
    
    % Shear from prestress (due to tendon inclination)
    % V = -P * de/dx (negative because upward force from draping)
    de_dx = gradient(e, x);
    V_prestress = V_prestress - Pe * de_dx;
end

% Normalize effective eccentricity by total force
total_P = sum(cellfun(@(t) t.Aps * t.fpi * (1 - prestress.losses), prestress.tendons));
if total_P > 0
    e_eff = e_eff / total_P;
end

end

%% INTERNAL FORCES - CORRECTED BOUNDARY CONDITIONS
function [V, M, V_dead, V_live, M_dead, M_live, reactions] = calculateInternalForces(beam, loads, section)
% Calculate shear and moment diagrams from applied loads
% CORRECTED: Proper boundary conditions for simple support

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
        w_dead(mask) = w_dead(mask) + w_segment;
    end
end

w = w_dead + w_live;

%% Calculate reactions from equilibrium
switch loads.support_type
    case 'simple'
        x_left = loads.supports(1);
        x_right = loads.supports(2);
        span = x_right - x_left;
        
        total_load = trapz(x, w);
        moment_left = trapz(x, w .* (x - x_left));
        
        for i = 1:size(loads.point, 1)
            x_p = loads.point(i, 1);
            P_p = loads.point(i, 2);
            M_p = loads.point(i, 3);
            
            total_load = total_load + P_p;
            moment_left = moment_left + P_p * (x_p - x_left) + M_p;
        end
        
        R_right = -moment_left / span;
        R_left = -total_load - R_right;
        
        reactions.left = R_left;
        reactions.right = R_right;
        reactions.x_left = x_left;
        reactions.x_right = x_right;
        
        fprintf('  Calculated reactions from equilibrium:\n');
        fprintf('    R_left  = %+.4f kips at x = %.1f in\n', R_left, x_left);
        fprintf('    R_right = %+.4f kips at x = %.1f in\n', R_right, x_right);
        fprintf('    Check: R_left + R_right + W_total = %.6f (should be ~0)\n', R_left + R_right + total_load);
        
    case 'cantilever'
        x_support = loads.supports(1);
        
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

%% Calculate shear and moment with CORRECT boundary conditions
switch loads.support_type
    case 'simple'
        [~, idx_left] = min(abs(x - reactions.x_left));
        [~, idx_right] = min(abs(x - reactions.x_right));
        
        V(idx_left) = reactions.left;
        M(idx_left) = 0;
        
        for i = (idx_left+1):n
            P_point = 0;
            M_point = 0;
            for j = 1:size(loads.point, 1)
                if abs(x(i) - loads.point(j, 1)) < dx/2
                    P_point = P_point + loads.point(j, 2);
                    M_point = M_point + loads.point(j, 3);
                end
            end
            
            w_avg = (w(i-1) + w(i)) / 2;
            V(i) = V(i-1) + w_avg * dx + P_point;
            
            V_avg = (V(i-1) + V(i)) / 2;
            M(i) = M(i-1) + V_avg * dx + M_point;
        end
        
        for i = (idx_left-1):-1:1
            P_point = 0;
            M_point = 0;
            for j = 1:size(loads.point, 1)
                if abs(x(i+1) - loads.point(j, 1)) < dx/2
                    P_point = P_point + loads.point(j, 2);
                    M_point = M_point + loads.point(j, 3);
                end
            end
            
            w_avg = (w(i) + w(i+1)) / 2;
            V(i) = V(i+1) - w_avg * dx - P_point;
            
            V_avg = (V(i) + V(i+1)) / 2;
            M(i) = M(i+1) - V_avg * dx - M_point;
        end
        
        fprintf('  Boundary condition verification:\n');
        fprintf('    M(0)   = %.6f kip-in (should be 0)\n', M(idx_left));
        fprintf('    M(L)   = %.6f kip-in (should be 0)\n', M(idx_right));
        fprintf('    V(0+)  = %+.4f kips (should equal R_left = %+.4f)\n', V(idx_left), reactions.left);
        fprintf('    V(L-)  = %+.4f kips (should equal -R_right = %+.4f)\n', V(idx_right), -reactions.right);
        
        if abs(M(idx_right)) > 0.1
            warning('M(L) = %.3f is not close to zero. Check reaction calculations!', M(idx_right));
        end
        
    case 'cantilever'
        [~, idx_support] = min(abs(x - reactions.x_support));
        
        V(idx_support) = reactions.vertical;
        M(idx_support) = reactions.moment;
        
        for i = idx_support+1:n
            P_point = 0;
            M_point = 0;
            for j = 1:size(loads.point, 1)
                if abs(x(i) - loads.point(j, 1)) < dx/2
                    P_point = P_point + loads.point(j, 2);
                    M_point = M_point + loads.point(j, 3);
                end
            end
            
            w_avg = (w(i-1) + w(i)) / 2;
            V(i) = V(i-1) + w_avg * dx + P_point;
            
            V_avg = (V(i-1) + V(i)) / 2;
            M(i) = M(i-1) + V_avg * dx + M_point;
        end
        
        for i = idx_support-1:-1:1
            P_point = 0;
            M_point = 0;
            for j = 1:size(loads.point, 1)
                if abs(x(i+1) - loads.point(j, 1)) < dx/2
                    P_point = P_point + loads.point(j, 2);
                    M_point = M_point + loads.point(j, 3);
                end
            end
            
            w_avg = (w(i) + w(i+1)) / 2;
            V(i) = V(i+1) - w_avg * dx - P_point;
            
            V_avg = (V(i) + V(i+1)) / 2;
            M(i) = M(i+1) - V_avg * dx - M_point;
        end
end

% Separate dead and live load effects (simplified approach)
V_dead = V;
M_dead = M;
V_live = zeros(1, n);
M_live = zeros(1, n);

end

%% STRESS CALCULATIONS
function stresses = calculateStresses(results, section, materials)
% Calculate stresses at critical locations
%
% Sign Convention:
%   Compression = negative
%   Tension = positive
%
% Prestress stress (using direct P, e approach — consistent with M_prestress):
%   f = -P/A + P*e*y/I   (where y is measured from centroid, + upward)
%   At top:    f_top = -P/A + P*e*yt/Ix    
%   At bottom: f_bot = -P/A - P*e*yb/Ix    
%
%   When e > 0 (tendon below): bottom gets MORE compression ✓
%
% External load stress:
%   f = -M*y/I  (for sagging M > 0: compression at top, tension at bottom)
%   At top:    f_top = -M*yt/Ix
%   At bottom: f_bot = +M*yb/Ix

n = length(results.x);

stresses.f_top = zeros(1, n);
stresses.f_bot = zeros(1, n);
stresses.f_top_prestress = zeros(1, n);
stresses.f_bot_prestress = zeros(1, n);
stresses.f_top_total = zeros(1, n);
stresses.f_bot_total = zeros(1, n);

A = section.A;
Ix = section.Ix;
yt = section.yt;
yb = section.yb;

for i = 1:n
    P = results.P(i);
    e = results.e(i);
    M = results.M(i);
    
    % Prestress stresses
    f_prestress_top = -P/A + P*e*yt/Ix;
    f_prestress_bot = -P/A - P*e*yb/Ix;
    
    % External load stresses
    f_load_top = -M*yt/Ix;
    f_load_bot = M*yb/Ix;
    
    stresses.f_top_prestress(i) = f_prestress_top;
    stresses.f_bot_prestress(i) = f_prestress_bot;
    stresses.f_top(i) = f_load_top;
    stresses.f_bot(i) = f_load_bot;
    stresses.f_top_total(i) = f_prestress_top + f_load_top;
    stresses.f_bot_total(i) = f_prestress_bot + f_load_bot;
end

stresses.fc_allow_compression = -0.45 * materials.fc;
stresses.fc_allow_tension = materials.fr;

end

%% CAPACITY CALCULATIONS
function capacity = calculateCapacity(section, materials, prestress, reinforcement)

fc = materials.fc;
fy = materials.fy;
fpu = materials.fpu;
fpy = materials.fpy;
beta1 = 0.85 - 0.05 * (fc - 4);
beta1 = max(0.65, min(0.85, beta1));

Aps_total = sum(cellfun(@(t) t.Aps, prestress.tendons));

yc = section.yc;
e_avg = mean(cellfun(@(t) mean(t.e), prestress.tendons));
dp = yc + e_avg;

fpe = prestress.tendons{1}.fpi * (1 - prestress.losses);

fps = fpu * (1 - 0.28 * (fpu / fc) * (Aps_total / (section.A * dp / section.yb)));
fps = min(fps, fpy);

fps_approx = min(fpe + 15, fpu);

b = section.A / (section.yt + section.yb);

c = Aps_total * fps_approx / (0.85 * fc * beta1 * b);
a = beta1 * c;

Mn = Aps_total * fps_approx * (dp - a/2);

epsilon_t = 0.003 * (dp - c) / c;
if epsilon_t >= 0.005
    phi = 0.9;
elseif epsilon_t <= 0.002
    phi = 0.65;
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
