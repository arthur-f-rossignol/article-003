%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%
%%     RESOURCE COMPETITION BETWEEN BUOYANCY-REGULATING AND SINKING      %%
%%         PHYTOPLANKTON SPECIES ALONG A STRATIFIED WATER COLUMN         %%
%%                                                                       %%
%%                 Arthur F. Rossignol & Sabine Wollrab                  %%
%%                                                                       %%
%%                               June 2026                               %%
%%                                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd "path/to/working/directory"

clearvars;
close all;

%% PARAMETERIZATION =======================================================

params = struct();

% Environmental parameters ------------------------------------------------

params.env.l   = 40;      % maximum depth of water column             [m]
params.env.l_t = 10;      % depth of thermocline                      [m]
params.env.l_b = 39;      % depth where the benthic layer begins      [m]
params.env.D_e = 100;     % eddy diffusion coefficient of epilimnion  [m²·day⁻¹]
params.env.D_h = 1;       % eddy diffusion coefficient of hypolimnion [m²·day⁻¹]
params.env.a_0 = 0.2;     % background turbidity                      [m⁻¹]
params.env.r   = 0.5;     % nutrient recycling rate
params.env.E   = 0.05;    % sediment interface permeability           [m⁻¹]
params.env.w_t = 0.5;     % width of thermocline                      [m]
params.env.N_0 = 50;      % nutrient concentration in the sediment    [mg(P)·m⁻³]
params.env.I_0 = 500;     % light intensity at surface                [μmol(photons)·m⁻²·s⁻¹]

% Buoyancy-regulating species' (A_BR) parameters --------------------------

params.sp1.mu  = 1;       % maximum growth rate                       [day⁻¹]
params.sp1.K   = 0.5;     % half-saturation constant for nutrient     [mg(P)·m⁻³]
params.sp1.H   = 120;     % half-saturation constant for light        [μmol(photons)·m⁻²·s⁻¹]
params.sp1.m   = 0.1;     % mortality rate                            [day⁻¹]
params.sp1.v   = 0.5;     % maximum vertical velocity                 [m·day⁻¹]
params.sp1.q   = 0.008;   % algal nutrient quota                      [mg(P)·mg(C)⁻¹]
params.sp1.a   = 0.0003;  % algal absorption coefficient              [m²·mg(C)⁻¹]
params.sp1.A_0 = 100;     % initial algal biomass                     [mg(C)·m⁻³]

% Sinking species' (A_S) parameters ---------------------------------------

params.sp2.mu  = 1;       % maximum growth rate                       [day⁻¹]
params.sp2.K   = 3;       % half-saturation constant for nutrient     [mg(P)·m⁻³]
params.sp2.H   = 20;      % half-saturation constant for light        [μmol(photons)·m⁻²·s⁻¹]
params.sp2.m   = 0.1;     % mortality rate                            [day⁻¹]
params.sp2.v   = 0.5;     % maximum vertical velocity                 [m·day⁻¹]
params.sp2.q   = 0.008;   % algal nutrient quota                      [mg(P)·mg(C)⁻¹]
params.sp2.a   = 0.0003;  % algal absorption coefficient              [m²·mg(C)⁻¹]
params.sp2.A_0 = 100;     % initial algal biomass                     [mg(C)·m⁻³]

% Numerical parameters ----------------------------------------------------

params.num.dz    = 0.08;  % spatial discretization step               [m]
params.num.t_max = 1e3;   % maximum simulation time                   [days]

%% SETUP ==================================================================

% Spatial discretization --------------------------------------------------

dz  = params.num.dz;
l   = params.env.l;
l_t = params.env.l_t;
l_b = params.env.l_b;

n   = floor(l / dz);
n_t = floor(l_t / dz);
Z   = linspace(0, l, n);

% Parameter vector --------------------------------------------------------

p = [dz; 
     n;
     n_t; 
     params.env.D_e; 
     params.env.D_h; 
     params.env.a_0; 
     params.env.r; 
     params.env.E; 
     params.env.w_t; 
     params.env.N_0; 
     params.env.I_0;
     params.sp1.mu;
     params.sp1.K; 
     params.sp1.H; 
     params.sp1.m; 
     params.sp1.v; 
     params.sp1.q; 
     params.sp1.a;
     params.sp2.mu; 
     params.sp2.K; 
     params.sp2.H; 
     params.sp2.m; 
     params.sp2.v; 
     params.sp2.q; 
     params.sp2.a];

% Initial conditions ------------------------------------------------------

U0 = [params.sp1.A_0 * ones(n, 1);                        
      params.sp2.A_0 * ones(n, 1);     
      zeros(n_t, 1); transpose(linspace(0, params.env.N_0, n - n_t))]; 

%% ODE SOLVING ============================================================

options = odeset('nonnegative', 1:(3 * n), 'RelTol', 1e-8, 'AbsTol', 1e-10);

fprintf('Starting simulation (t_max = %.0f days)...\n', params.num.t_max);

[t, sol] = ode15s(@(t, U) equations(t, U, p), [0, params.num.t_max], U0, options);

fprintf('Simulation completed.\n');

%% POST-PROCESSING ========================================================

% Equilibrium profiles ----------------------------------------------------

A1 = sol(end, 1:n);
A2 = sol(end, (n + 1):(2 * n));
N  = sol(end, (2 * n + 1):(3 * n));

% Equilibrium biomass -----------------------------------------------------

biomass = compute_biomass_distribution(A1, A2, dz, l_t, l_b, l, n);

fprintf('\nBiomass distribution summary:');
fprintf('\n  Buoyancy-regulating species (A_BR):\n');
fprintf('    epilimnion:    %.2f\n', biomass.sp1.epilimnion);
fprintf('    hypolimnion:   %.2f\n', biomass.sp1.hypolimnion);
fprintf('    benthic layer: %.2f\n', biomass.sp1.benthic);
fprintf('    total:         %.2f\n', biomass.sp1.total);
fprintf('\n  Sinking species (A_S):\n');
fprintf('    epilimnion:    %.2f\n', biomass.sp2.epilimnion);
fprintf('    hypolimnion:   %.2f\n', biomass.sp2.hypolimnion);
fprintf('    benthic layer: %.2f\n', biomass.sp2.benthic);
fprintf('    total:         %.2f\n', biomass.sp2.total);

% Plotting ----------------------------------------------------------------

profiles = compute_profiles(A1, A2, N, p);
plot_equilibrium_profiles(A1, A2, N, profiles, Z, dz, n);
plot_temporal_evolution(t, sol, Z, n);

%% HELPER FUNCTIONS =======================================================

% Computing biomass distribution ------------------------------------------

function biomass = compute_biomass_distribution(A1, A2, dz, l_t, l_b, l, n)
    
    % boundary indices
    idx_thermocline = max(2, min(round(l_t / dz) + 1,  n + 1));
    idx_benthic     = max(idx_thermocline, min(round(l_b / dz) + 1, n + 1));

    % epilimnion
    ep = 1:(idx_thermocline - 1);
    biomass.sp1.epilimnion = sum(A1(ep)) * dz;
    biomass.sp2.epilimnion = sum(A2(ep)) * dz;

    % hypolimnion
    hy = idx_thermocline:(idx_benthic - 1);
    biomass.sp1.hypolimnion = sum(A1(hy)) * dz;
    biomass.sp2.hypolimnion = sum(A2(hy)) * dz;

    % benthic layer
    be = idx_benthic:n;
    biomass.sp1.benthic = sum(A1(be)) * dz;
    biomass.sp2.benthic = sum(A2(be)) * dz;

    % total biomass
    biomass.sp1.total = biomass.sp1.epilimnion ...
                      + biomass.sp1.hypolimnion ...
                      + biomass.sp1.benthic;
    biomass.sp2.total = biomass.sp2.epilimnion ...
                      + biomass.sp2.hypolimnion ...
                      + biomass.sp2.benthic;
end

% Computing profiles ------------------------------------------------------

function profiles = compute_profiles(A1, A2, N, p)
    
    % parameter unpacking
    dz   = p(1);
    n    = p(2);
    n_t  = p(3);
    D_e  = p(4);
    D_h  = p(5);
    a_0  = p(6);
    w_t  = p(9);
    I_0  = p(11);
    mu_1 = p(12);
    K_1  = p(13);
    H_1  = p(14);
    m_1  = p(15);
    v_1  = p(16);
    a_1  = p(18);
    mu_2 = p(19);
    K_2  = p(20);
    H_2  = p(21);
    m_2  = p(22);
    v_2  = p(23);
    a_2  = p(25);
    
    % preallocation of vectors
    profiles.I      = zeros(n, 1);
    profiles.dG1_dz = zeros(n, 1);
    profiles.dG2_dz = zeros(n, 1);
    
    % light intensity
    profiles.I(1) = I_0 * exp(- a_0 * 0.5 * dz ...
                              - a_1 * ((3 * A1(1) - A1(2)) / 8 + 3 * A1(1) / 4) * dz ...
                              - a_2 * ((3 * A2(1) - A2(2)) / 8 + 3 * A2(1) / 4) * dz);
    for i = 2:n
        M = - a_1 * ((3 * A1(1) - A1(2)) / 8 + 3 * A1(1) / 4) * dz ...
            - a_2 * ((3 * A2(1) - A2(2)) / 8 + 3 * A2(1) / 4) * dz;
        for k = 2:(i - 1)
            M = M - a_1 * A1(k) * dz ...
                  - a_2 * A2(k) * dz;
        end
        profiles.I(i) = I_0 * exp(M - a_0 * (i - 0.5) * dz ...
                                    - a_1 * A1(i) * 0.5 * dz ...
                                    - a_2 * A2(i) * 0.5 * dz);
    end
    
    % growth rates
    profiles.G1 = mu_1 * min(N(:) ./ (K_1 + N(:)), profiles.I ./ ...
                  (H_1 + profiles.I)) - m_1;
    profiles.G2 = mu_2 * min(N(:) ./ (K_2 + N(:)), profiles.I ./ ...
                  (H_2 + profiles.I)) - m_2;
    
    % fitness gradients
    profiles.dG1_dz(2:(n - 1)) = (profiles.G1(3:n) - profiles.G1(1:(n - 2))) / (2 * dz);
    profiles.dG1_dz(1)         = (profiles.G1(2) - profiles.G1(1)) / dz;
    profiles.dG1_dz(n)         = (profiles.G1(n) - profiles.G1(n - 1)) / dz;
    profiles.dG2_dz(2:(n - 1)) = (profiles.G2(3:n) - profiles.G2(1:(n - 2))) / (2 * dz);
    profiles.dG2_dz(1)         = (profiles.G2(2) - profiles.G2(1)) / dz;
    profiles.dG2_dz(n)         = (profiles.G2(n) - profiles.G2(n - 1)) / dz;

    % vertical velocities
    profiles.V1 = v_1 * profiles.dG1_dz ./ (abs(profiles.dG1_dz) + 1e-3);
    profiles.V2 = v_2 * ones(n, 1);

    % eddy diffusion
    profiles.D = D_h + (D_e - D_h) ./ (1 + exp(((1:n)' - n_t) * dz / w_t));
end

% Plotting equilibrium profiles -------------------------------------------

function plot_equilibrium_profiles(A1, A2, N, profiles, Z, dz, n)
    
fig = figure('Position', [100, 100, 1400, 900]);
    
    depth = linspace(0, n * dz, n);
    
    % A_BR - biomass
    subplot(3, 4, 1);
    plot(A1, depth, 'LineWidth', 2.5);
    title('algal biomass of A_{BR} [mg(C)·m⁻³]');
    ylabel('depth [m]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % A_BR - growth
    subplot(3, 4, 2);
    plot(profiles.G1, depth, 'LineWidth', 2.5);
    title('net growth rate of A_{BR} [m⁻¹]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % A_BR - fitness gradient
    subplot(3, 4, 3);
    plot(profiles.dG1_dz, depth, 'LineWidth', 2.5);
    title('fitness gradient of A_{BR} [m⁻²]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % A_BR - velocity
    subplot(3, 4, 4);
    plot(profiles.V1, depth, 'LineWidth', 2.5);
    title('vertical velocity of A_{BR} [m·day⁻¹]');
    set(gca, 'YDir', 'reverse');
    grid on;
    
    % A_S - biomass
    subplot(3, 4, 5);
    plot(A2, depth, 'LineWidth', 2.5);
    title('algal biomass of A_S [mg(C)·m⁻³]');
    ylabel('depth [m]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % A_S - growth
    subplot(3, 4, 6);
    plot(profiles.G2, depth, 'LineWidth', 2.5);
    title('net growth rate of A_S [m⁻¹]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % A_S - fitness gradient
    subplot(3, 4, 7);
    plot(profiles.dG2_dz, depth, 'LineWidth', 2.5);
    title('fitness gradient of A_S [m⁻²]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % A_S - velocity
    subplot(3, 4, 8);
    plot(profiles.V2, depth, 'LineWidth', 2.5);
    title('vertical velocity of A_S [m·day⁻¹]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % total biomass
    subplot(3, 4, 9);
    plot(A1 + A2, depth, 'LineWidth', 2.5);
    title('total algal biomass [mg(C)·m⁻³]');
    ylabel('Depth [m]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % nutrients
    subplot(3, 4, 10);
    plot(N, depth, 'LineWidth', 2.5);
    title('nutrient concentration [mg(P)·m⁻³]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % light
    subplot(3, 4, 11);
    plot(profiles.I, depth, 'LineWidth', 2.5);
    title('light intensity [μmol(photons)·m⁻²·s⁻¹]');
    set(gca, 'YDir','reverse');
    grid on;
    
    % diffusion
    subplot(3, 4, 12);
    plot(profiles.D, depth, 'LineWidth', 2.5);
    title('eddy diffusion coefficient [m²·day⁻¹]');
    set(gca, 'YDir','reverse');
    grid on;
    
    saveas(fig, 'equilibrium.fig');
end

% Plotting temporal evolution ---------------------------------------------

function plot_temporal_evolution(t, sol, Z, n)

    figure('Position', [100, 100, 1200, 500]);
    
    subplot(1, 2, 1);
    contourf(t, Z, transpose(sol(:, 1:n)), 80, 'LineColor', 'none');
    colorbar;
    title('buoyancy-regulating species - temporal evolution');
    xlabel('Time [days]');
    ylabel('Depth [m]');
    set(gca, 'YDir', 'reverse');
    
    subplot(1, 2, 2);
    contourf(t, Z, transpose(sol(:, (n + 1):(2 * n))), 80, 'LineColor', 'none');
    colorbar;
    title('sinking species - temporal evolution');
    xlabel('Time [days]');
    ylabel('Depth [m]');
    set(gca, 'YDir', 'reverse');
end

% Parameter sweep (example: hypolimnetic eddy diffusion) ------------------

function results = parameter_sweep_diffusion(params, U0, options)
    
    step = 100;
    PARAM = logspace(log10(0.05), log10(50), step);
    
    dz  = params.num.dz;
    l   = params.env.l;
    l_t = params.env.l_t;
    n   = floor(l / dz);
    n_t = floor(l_t / dz);
    
    options = odeset('nonnegative', 1:(3 * n), 'RelTol', 1e-8, 'AbsTol', 1e-10);
    
    results.PARAM = PARAM;
    results.Z     = linspace(0, l, n);
    results.A1    = NaN(n, length(PARAM));
    results.A2    = NaN(n, length(PARAM));
    results.N     = NaN(n, length(PARAM));
    
    for i = 1:length(PARAM)
    
        D_h = PARAM(i);
        
        p = [dz; 
             n; 
             n_t; 
             params.env.D_e; 
             D_h; 
             params.env.a_0; 
             params.env.r; 
             params.env.E; 
             params.env.w_t; 
             params.env.N_0; 
             params.env.I_0;
             params.sp1.mu; 
             params.sp1.K; 
             params.sp1.H; 
             params.sp1.m; 
             params.sp1.v; 
             params.sp1.q; 
             params.sp1.a;
             params.sp2.mu; 
             params.sp2.K; 
             params.sp2.H; 
             params.sp2.m; 
             params.sp2.v; 
             params.sp2.q; 
             params.sp2.a];
        
        [~, sol] = ode15s(@(t, U) equations(t, U, p), [0, params.num.t_max], U0, options);
        
        results.A1(:, i) = sol(end, 1:n);
        results.A2(:, i) = sol(end, (n + 1):(2 * n));
        results.N(:, i)  = sol(end, (2 * n + 1):(3 * n));
    end
end

% Plot parameter sweep ----------------------------------------------------

function plot_parameter_sweep(results, l)
   
    fig = figure('Position', [100, 100, 1200, 500]);
    
    MAX = max(max(max(results.A1)), max(max(results.A2)));
    MIN = min(min(min(results.A1)), min(min(results.A2)));
    
    subplot(1, 2, 1);
    contourf(results.PARAM, fliplr(results.Z), flipud(results.A1), 30, 'LineColor', 'none');
    caxis([MIN MAX]);
    colorbar;
    title('A_BR [mg(C)·m⁻³]');
    ylabel('depth [m]');
    xlabel('eddy diffusion coefficient [m²·day⁻¹]');
    set(gca, 'YDir', 'reverse', 'XScale', 'log');
    
    subplot(1, 2, 2);
    contourf(results.PARAM, fliplr(results.Z), flipud(results.A2), 30, 'LineColor', 'none');
    caxis([MIN MAX]);
    colorbar;
    title('A_S [mg(C)·m⁻³]');
    ylabel('depth [m]');
    xlabel('eddy diffusion coefficient [m²·day⁻¹]');
    set(gca, 'YDir', 'reverse', 'XScale', 'log');
    
    saveas(fig, 'parameter_sweep.fig');
end

% Equations ---------------------------------------------------------------

function dU_dt = equations(t, U, p)
    
    % parameter unpacking
    dz   = p(1);
    n    = p(2);
    n_t  = p(3);
    D_e  = p(4);
    D_h  = p(5);
    a_0  = p(6);
    r    = p(7);
    E    = p(8);
    w_t  = p(9);
    N_0  = p(10);
    I_0  = p(11);
    mu_1 = p(12);
    K_1  = p(13);
    H_1  = p(14);
    m_1  = p(15);
    v_1  = p(16);
    q_1  = p(17);
    a_1  = p(18);
    mu_2 = p(19);
    K_2  = p(20);
    H_2  = p(21);
    m_2  = p(22);
    v_2  = p(23);
    q_2  = p(24);
    a_2  = p(25);

    % preallocation of vectors
    I        = zeros(n, 1);
    dA1_dt   = zeros(n, 1);
    dV1A1_dz = zeros(n, 1);
    dDA1_dz2 = zeros(n, 1);
    dG1_dz   = zeros(n, 1);
    dA2_dz   = zeros(n, 1);
    dA2_dt   = zeros(n, 1);
    dDA2_dz2 = zeros(n, 1);
    dN_dt    = zeros(n, 1);
    dDN_dz2  = zeros(n, 1);

    % starting values
    A1 = U(1:n);
    A2 = U((n + 1):(2 * n));
    N  = U((2 * n + 1):(3 * n));

    % diffusion
    D  = D_h + (D_e - D_h) ./ (1 + exp(((1:n)' - n_t) * dz / w_t));
    Df = 0.5 * (D(1:(n - 1)) + D(2:n));

    % light intensity
    I(1) = I_0 * exp(- a_0 * 0.5 * dz ...
                     - a_1 * ((3 * A1(1) - A1(2)) / 8 + 3 * A1(1) / 4) * dz ...
                     - a_2 * ((3 * A2(1) - A2(2)) / 8 + 3 * A2(1) / 4) * dz);
    for i = 2:n
        M = - a_1 * ((3 * A1(1) - A1(2)) / 8 + 3 * A1(1) / 4) * dz ...
            - a_2 * ((3 * A2(1) - A2(2)) / 8 + 3 * A2(1) / 4) * dz;
        for k = 2:(i - 1)
            M = M - a_1 * A1(k) * dz ...
                  - a_2 * A2(k) * dz;
        end
        I(i) = I_0 * exp(M - a_0 * (i - 0.5) * dz ...
                           - a_1 * A1(i) * 0.5 * dz ...
                           - a_2 * A2(i) * 0.5 * dz);
    end

    % growth rates
    G1 = mu_1 * min(N ./ (K_1 + N), I ./ (H_1 + I));
    G2 = mu_2 * min(N ./ (K_2 + N), I ./ (H_2 + I));
  
    % fitness gradient of A_BR 
    dG1_dz(2:(n - 1)) = (G1(3:n) - G1(1:(n - 2))) / (2 * dz);
    dG1_dz(1)         = (G1(2) - G1(1)) / dz;
    dG1_dz(n)         = (G1(n) - G1(n - 1)) / dz;

    % vertical velocity of A_BR
    V1  = v_1 * dG1_dz ./ (abs(dG1_dz) + 1e-3);
    V1f = 0.5 * (V1(1:(n - 1)) + V1(2:n));

    % upwind second-order scheme for advection term of A_BR's biomass 
    ii = 3:(n - 2);
    V1f_r = V1f(ii); 
    V1f_l = V1f(ii - 1); 
    tr = (V1f_r > 0);
    tl = (V1f_l > 0);
    dV1A1_dz(ii) = (V1f_r .* (2 * A1(ii) + 5 * A1(ii + 1) - A1(ii + 2)) .* (1 - tr) ...
                 + V1f_r .* (- A1(ii - 1) + 5 * A1(ii) + 2 * A1(ii + 1)) .* tr ...
                 - V1f_l .* ( 2 * A1(ii - 1) + 5 * A1(ii) - A1(ii + 1)) .* (1 - tl) ...
                 - V1f_l .* (- A1(ii - 2) + 5 * A1(ii - 1) + 2 * A1(ii)) .* tl) / (6 * dz);
    t1 = (V1f(1) > 0);
    dV1A1_dz(1) = V1f(1) * t1 * (A1(1) + A1(2)) / (2 * dz) ...
                + V1f(1) * (1 - t1) * (2 * A1(1) + 5 * A1(2) - A1(3)) / (6 * dz);
    t1 = (V1f(2) > 0);
    t2 = (V1f(1) > 0);
    dV1A1_dz(2) = ((- V1f(2) * t1 - 3 * V1f(1) * t2 - 2 * V1f(1) * (1 - t2)) * A1(1) ...
                + (2 * V1f(2) * (1 - t1) + 5 * V1f(2) * t1 ...
                - 5 * V1f(1) * (1 - t2) - 3 * V1f(1) * t2) * A1(2) ...
                + (5 * V1f(2) * (1 - t1) + 2 * V1f(2) * t1 + V1f(1) * (1 - t2)) * A1(3) ...
                - V1f(2) * (1 - t1) * A1(4)) / (6 * dz);
    t1 = (V1f(n - 1) > 0);
    t2 = (V1f(n - 2) > 0);
    dV1A1_dz(n - 1) = (V1f(n - 2) * t2 * A1(n - 3) ...
                    + (- V1f(n - 1) * t1 - 2 * V1f(n - 2) * (1 - t2) ...
                    - 5 * V1f(n - 2) * t2) * A1(n - 2) ...
                    + (5 * V1f(n - 1) * t1 + 3 * V1f(n - 1) * (1 - t1) ...
                    - 5 * V1f(n - 2) * (1 - t2) - 2 * V1f(n - 2) * t2) * A1(n - 1) ...
                    + (2 * V1f(n - 1) * t1 + 3 * V1f(n - 1) * (1 - t1) ...
                    + V1f(n - 2) * (1 - t2)) * A1(n)) / (6 * dz);
    t1 = (V1(n) > 0);
    t2 = (V1f(n - 1) > 0);
    dV1A1_dz(n) = V1(n) * t1 * (- A1(n - 1) + 7 * A1(n)) / (6 * dz) ...
                + V1(n) * (1 - t1) * A1(n) / dz ...
                - V1f(n - 1) * (1 - t2) * (A1(n) + A1(n - 1)) / (2 * dz) ...
                - V1f(n - 1) * t2 * (2 * A1(n) + 5 * A1(n - 1) - A1(n - 2)) / (6 * dz);

    % upwind second-order scheme for advection term of A_S's biomass 
    dA2_dz(3:(n - 1)) = (2 * A2(4:n) ...
                      + 3 * A2(3:(n - 1)) ...
                      - 6 * A2(2:(n - 2)) ...
                      + A2(1:(n - 3))) / (6 * dz); 
    dA2_dz(1)         = (A2(1) + A2(2)) / (2 * dz);
    dA2_dz(2)         = (- A2(1) + 5 * A2(2) + 2 * A2(3)) / (6 * dz) ...
                      - (A2(1) + A2(2)) / (2 * dz);
    dA2_dz(n)         = (A2(n - 2) - 6 * A2(n - 1) + 5 * A2(n)) / (6 * dz);
    
    % symmetric 2nd-order scheme for diffusion term of A_BR's biomass
    dDA1_dz2(2:(n - 1)) = (Df(2:(n - 1)) .* (A1(3:n)   - A1(2:(n - 1))) ...
                        - Df(1:(n - 2)) .* (A1(2:(n - 1)) - A1(1:(n - 2))) ) / dz^2;
    dDA1_dz2(1)         = Df(1) * (A1(2) - A1(1)) / dz^2;
    dDA1_dz2(n)         = - Df(n - 1) * (A1(n) - A1(n - 1)) / dz^2;
    
    % symmetric 2nd-order scheme for diffusion term of A_S's biomass
    dDA2_dz2(2:(n - 1)) = (Df(2:(n - 1)) .* (A2(3:n) - A2(2:(n - 1))) ...
                        - Df(1:(n - 2)) .* (A2(2:(n - 1)) - A2(1:(n - 2)))) / dz^2;
    dDA2_dz2(1)         = Df(1) * (A2(2) - A2(1)) / dz^2;
    dDA2_dz2(n)         = - Df(n - 1) * (A2(n) - A2(n - 1)) / dz^2;

    % symmetric 2nd-order scheme for diffusion term of nutrients
    dDN_dz2(2:(n - 1)) = (Df(2:(n - 1)) .* (N(3:n) - N(2:(n - 1))) ...
                       - Df(1:(n - 2)) .* (N(2:(n - 1)) - N(1:(n - 2))) ) / dz^2;
    dDN_dz2(1)         = Df(1) * (N(2) - N(1)) / dz^2;
    dDN_dz2(n)         = D(n) * E * (N_0 - N(n)) / dz ...
                       - Df(n - 1) * (N(n) - N(n - 1)) / dz^2;
                       
    % time derivatives of A_BR's biomass, A_S's biomass, nutrients
    dA1_dt(1:n) = (G1(1:n) - m_1) .* A1(1:n) - dV1A1_dz(1:n) + dDA1_dz2(1:n);
    dA2_dt(1:n) = (G2(1:n) - m_2) .* A2(1:n) - v_2 * dA2_dz(1:n) + dDA2_dz2(1:n);
    dN_dt(1:n)  = q_1 * (r * m_1 - G1(1:n)) .* A1(1:n) ...
                + q_2 * (r * m_2 - G2(1:n)) .* A2(1:n) ...
                + dDN_dz2(1:n);
    
    % output
    dU_dt = [dA1_dt; dA2_dt; dN_dt];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
