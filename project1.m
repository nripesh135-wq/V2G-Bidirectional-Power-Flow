%% 
%% MPC_V2G_Controller_runready.m
% Run-ready MPC for multi-EV bidirectional power flow (copy-paste into MATLAB)
% Requires Optimization Toolbox (quadprog)
% Output: SOC plots, per-EV power, net grid load, cost
clear; close all; clc;

%% ---------------- PARAMETERS ----------------
Tsim = 24;       % total simulation steps (e.g., hours)
dt = 1;          % hour per step
N = 6;           % MPC horizon (steps)
nEV = 3;         % number of EVs

% EV parameters (column vectors)
Ebat = [40; 60; 30];        % kWh
Pmax_ch = [7; 7; 3.3];      % kW
Pmax_dis = [7; 7; 3.3];     % kW
eta_ch = 0.95;
eta_dis = 0.95;
SOC_min = 0.2; SOC_max = 0.9;
SOC0_frac = [0.5; 0.8; 0.3];  % initial SOC fraction

% Objective weights (tune these)
w_grid = 20;    % weight for grid smoothing
w_cost = 0.5;   % charging cost weight (linear)
w_rev = 0.6;    % discharge revenue weight (linear, reward)
w_deg = 0.01;   % degradation quadratic weight per power^2
w_small = 1e-3; % small coupling to discourage simultaneous ch+dis

% Reference load (target)
L_ref = 300;    % kW

% Synthetic base load and price (demo). Replace with CSV if you want.
base_load = 300 + 50*sin((1:Tsim)/Tsim*2*pi);   % kW (1 x Tsim)
price = 0.10 + 0.07*(sin((1:Tsim)/Tsim*2*pi)+1); % $/kWh (1 x Tsim)
p_rev = price;  % revenue price

%% ---------------- PREALLOCATE LOGS ----------------
SOC = zeros(nEV, Tsim+1); SOC(:,1) = SOC0_frac .* Ebat; % energy in kWh
Pch = zeros(nEV, Tsim);
Pdis = zeros(nEV, Tsim);
grid_net = zeros(1, Tsim);
cost_log = zeros(1, Tsim);

options = optimoptions('quadprog','Display','off','TolFun',1e-9);

%% ---------------- MPC LOOP ----------------
for t = 1:Tsim
    % Horizon indices
    idx = t : min(t+N-1, Tsim);
    Nh = length(idx);
    d = base_load(idx)';   % Nh x 1
    c = price(idx)';       % Nh x 1
    r = p_rev(idx)';       % Nh x 1

    nv_per_step = 2 * nEV;
    Udim = Nh * nv_per_step;   % total decision vars

    % Initialize H and f for quadprog (quadprog minimizes 0.5*x'H*x + f'x)
    H = zeros(Udim, Udim);
    f = zeros(Udim, 1);

    % Build H and f
    for k = 1:Nh
        off = (k-1)*nv_per_step;
        idx_ch = off + (1:nEV);
        idx_dis = off + (nEV+1:2*nEV);

        % a vector mapping to net EV injection at this step: sum(ch) - sum(dis)
        a = zeros(1, Udim);
        a(idx_ch) = 1;
        a(idx_dis) = -1;

        % grid smoothing quadratic: w_grid * (a*U + (d-Lref))^2
        % Quadratic part contributes: 2*w_grid * (a' * a) (because quadprog has 1/2 factor)
        H = H + 2 * w_grid * (a' * a);

        % linear part from expansion: 2*w_grid*(d-Lref)*a' (added to f)
        c0 = d(k) - L_ref;
        f = f + 2 * w_grid * c0 * a';    % f is column

        % linear cost/revenue on charging/discharging (energy price * dt)
        f(idx_ch) = f(idx_ch) + w_cost * c(k) * dt;   % charging increases cost
        f(idx_dis) = f(idx_dis) - w_rev * r(k) * dt;  % discharging gives revenue (reduces objective)

        % degradation: w_deg * (p_ch^2 + p_dis^2)  => contributes 2*w_deg on diagonal (quadprog half factor)
        for j = 1:nEV
            H(idx_ch(j), idx_ch(j)) = H(idx_ch(j), idx_ch(j)) + 2 * w_deg;
            H(idx_dis(j), idx_dis(j)) = H(idx_dis(j), idx_dis(j)) + 2 * w_deg;
            % small coupling to discourage simultaneous ch + dis: we add off-diagonal sym term
            H(idx_ch(j), idx_dis(j)) = H(idx_ch(j), idx_dis(j)) + w_small;
            H(idx_dis(j), idx_ch(j)) = H(idx_dis(j), idx_ch(j)) + w_small;
        end
    end

    % numeric regularization to ensure positive definiteness
    H = H + 1e-6 * eye(Udim);

    % Bounds (lb and ub) for each variable
    lb = zeros(Udim,1);
    ub = inf(Udim,1);
    for k = 1:Nh
        off = (k-1)*nv_per_step;
        idx_ch = off + (1:nEV);
        idx_dis = off + (nEV+1:2*nEV);
        lb(idx_ch) = 0;
        ub(idx_ch) = Pmax_ch;
        lb(idx_dis) = 0;
        ub(idx_dis) = Pmax_dis;
    end

    % SOC constraints (linear inequalities): G * U <= h
    % Use energy form: E_{j, k+ell} = Ecur + dt * sum_{i=0..ell-1} (eta_ch*p_ch - (1/eta_dis)*p_dis)
    % For each EV j and each future step ell, enforce:
    %   Ecur + coeff*U <= SOC_max*Ebat
    %  -Ecur - coeff*U <= -SOC_min*Ebat   => (-coeff)*U <= Ecur - SOC_min*Ebat
    G = [];
    h = [];
    for j = 1:nEV
        Ecur = SOC(j, t);  % current stored energy (kWh)
        for ell = 1:Nh
            coeff = zeros(1, Udim);
            for kk = 1:ell
                off2 = (kk-1)*nv_per_step;
                coeff(off2 + j) = coeff(off2 + j) + eta_ch * dt;            % p_ch contribution (kWh per kW)
                coeff(off2 + nEV + j) = coeff(off2 + nEV + j) - (1/eta_dis) * dt; % p_dis contribution
            end
            % Upper bound: coeff * U <= SOC_max*Ebat(j) - Ecur
            G = [G; coeff];
            h = [h; SOC_max * Ebat(j) - Ecur];
            % Lower bound: -coeff * U <= Ecur - SOC_min*Ebat(j)
            G = [G; -coeff];
            h = [h; Ecur - SOC_min * Ebat(j)];
        end
    end

    % Solve QP: minimize 0.5 U'HU + f'U subject to G*U <= h and lb <= U <= ub
    % quadprog syntax: quadprog(H, f, A, b, Aeq, beq, lb, ub)
    try
        [Uopt, fval, exitflag] = quadprog(H, f, G, h, [], [], lb, ub, [], options);
    catch ME
        warning('quadprog error at t=%d: %s', t, ME.message);
        Uopt = zeros(Udim,1);
        exitflag = -1;
    end

    if isempty(Uopt)
        Uopt = zeros(Udim,1);
    end

    % Extract first-step controls
    Ustep = Uopt(1:nv_per_step);
    pch_step = Ustep(1:nEV);
    pdis_step = Ustep(nEV+1:2*nEV);

    % Apply and update SOC
    for j = 1:nEV
        deltaE = (eta_ch * pch_step(j) - (1/eta_dis) * pdis_step(j)) * dt; % kWh change
        SOC(j, t+1) = SOC(j, t) + deltaE;
        % Clip safety
        SOC(j, t+1) = min(max(SOC(j, t+1), SOC_min * Ebat(j)), SOC_max * Ebat(j));
        Pch(j, t) = pch_step(j);
        Pdis(j, t) = pdis_step(j);
    end

    % Net grid load after EVs
    grid_net(t) = d(1) + sum(pch_step) - sum(pdis_step);

    % Step cost (reporting)
    cost_log(t) = w_cost * price(t) * sum(pch_step) * dt - w_rev * p_rev(t) * sum(pdis_step) * dt ...
                  + w_deg * sum(pch_step.^2 + pdis_step.^2);
end

%% ---------------- PLOTS ----------------
SOC_frac = SOC ./ Ebat; % fraction

figure('Name','SOC Trajectories');
plot(0:Tsim, SOC_frac', 'LineWidth', 1.5); xlabel('Time step'); ylabel('SOC (fraction)');
legend(arrayfun(@(j) sprintf('EV%d', j), 1:nEV, 'UniformOutput', false)); grid on; title('SOC vs Time');

figure('Name','Per-EV Power');
plot(1:Tsim, Pch', '--', 1:Tsim, Pdis'); xlabel('Time step'); ylabel('Power (kW)');
legend([arrayfun(@(j) sprintf('EV%d ch', j), 1:nEV, 'UniformOutput', false), arrayfun(@(j) sprintf('EV%d dis', j), 1:nEV, 'UniformOutput', false)]);
grid on; title('Per-EV charge (dashed) and discharge (solid)');

figure('Name','Grid Net Load');
plot(1:Tsim, base_load, 'r--', 1:Tsim, grid_net, 'r-', 'LineWidth', 1.2); xlabel('Time step'); ylabel('kW');
legend('Base load', 'Net load after EVs'); grid on; title('Grid net load');

figure('Name','Price and Cost');
yyaxis left; plot(1:Tsim, price, '-o'); ylabel('Price ($/kWh)');
yyaxis right; plot(1:Tsim, cost_log, '-x'); ylabel('Step cost (units)');
xlabel('Time step'); legend('Price', 'Step cost'); grid on; title('Price and Step Cost');

fprintf('Simulation finished. Total reported cost (sum step costs): %f\n', sum(cost_log));
