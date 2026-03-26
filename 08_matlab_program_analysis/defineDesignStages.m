function stages = defineDesignStages(materials, prestress)
% DEFINEDESIGNSTAGES  Define ACI 318-19 design stages for prestressed beam analysis.
%
% USER TOGGLE: Set stage.active = false to skip that stage.
%
% Sign convention: compression = POSITIVE, tension = NEGATIVE (ACI 318 / Naaman)
%
% Inputs:
%   materials  - materials struct from inputPrestressedBeam_*() containing:
%                .f_ci_allow, .f_ti_allow, .f_cs_allow_sust,
%                .f_cs_allow_total, .f_ts_allow, .fc, .fci
%   prestress  - prestress struct containing .losses (e.g. 0.15 = 15%)
%
% Output:
%   stages  - cell array of stage structs, each with fields:
%     .name          string label used in titles and saved filenames
%     .eta           effective prestress factor (1 - fraction of losses realized)
%                      Transfer     : eta = 1.0  (no losses yet at release)
%                      Service      : eta = 1 - prestress.losses  (e.g. 0.85)
%     .include_sw    logical: include self-weight load
%     .include_sdl   logical: include superimposed dead load
%     .include_ll    logical: include live load
%     .f_allow_compr allowable compression stress [ksi, positive]
%     .f_allow_tens  allowable tension stress [ksi, negative]
%     .active        logical: compute this stage (set to false to skip)

eta_service = 1 - prestress.losses;   % e.g. 0.85 when losses = 0.15

%% ---- Stage 1: Transfer (immediately after prestress release) ----
%   - Prestress: full initial force (no losses yet)   → eta = 1.0
%   - Loads: self-weight only (beam on dunnage / being lifted)
%   - Allowable: ACI 318-19 Table 24.5.3.1 (using f'ci)
s1.name          = 'Transfer';
s1.eta           = 1.0;
s1.include_sw    = true;
s1.include_sdl   = false;
s1.include_ll    = false;
s1.f_allow_compr = materials.f_ci_allow;    % +0.60*f'ci  (e.g. +2.88 ksi)
s1.f_allow_tens  = materials.f_ti_allow;    % -3√f'ci_psi/1000 (e.g. -0.208 ksi)
s1.active        = true;   % ← USER TOGGLE: set false to skip Transfer

%% ---- Stage 2: Service — Sustained loads (SW + SDL, no LL) ----
%   - Prestress: after losses   → eta = 1 - losses
%   - Loads: self-weight + superimposed dead load
%   - Allowable: ACI 318-19 Table 24.5.4.1 sustained compression limit
s2.name          = 'Service_Sustained';
s2.eta           = eta_service;
s2.include_sw    = true;
s2.include_sdl   = true;
s2.include_ll    = false;
s2.f_allow_compr = materials.f_cs_allow_sust;   % +0.45*f'c  (e.g. +2.70 ksi)
s2.f_allow_tens  = materials.f_ts_allow;         % -12√f'c_psi/1000 Class C (e.g. -0.929 ksi)
s2.active        = true;   % ← USER TOGGLE: set false to skip Service_Sustained

%% ---- Stage 3: Service — Total loads (SW + SDL + LL) ----
%   - Prestress: after losses   → eta = 1 - losses
%   - Loads: all loads combined
%   - Allowable: ACI 318-19 Table 24.5.4.1 total load compression limit
s3.name          = 'Service_Total';
s3.eta           = eta_service;
s3.include_sw    = true;
s3.include_sdl   = true;
s3.include_ll    = true;
s3.f_allow_compr = materials.f_cs_allow_total;  % +0.60*f'c  (e.g. +3.60 ksi)
s3.f_allow_tens  = materials.f_ts_allow;
s3.active        = true;   % ← USER TOGGLE: set false to skip Service_Total

stages = {s1, s2, s3};

end
