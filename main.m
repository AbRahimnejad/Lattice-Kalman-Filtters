% simulate_eha_lkf.m
% ============================================================================
% Lattice Kalman Filter (LKF) demo on a nonlinear Electro-Hydrostatic Actuator
% (EHA) model (discrete-time) with full-state measurement for clarity.
%
% This script:
%   1) defines a stiff, mildly non-linear 4-state plant from the literature;
%   2) simulates the "true" system under process/measurement noise;
%   3) runs a Lattice Kalman Filter (LKF) based on quasi–Monte Carlo (Korobov)
%      sampling pushed through the Gaussian inverse CDF (normal score mapping);
%   4) reports RMSE and plots truth vs. estimate and estimation error.
%
% Files required in the path:
%   - lkf.m  (contains the LKF implementation and helpers)
%
% Notes on numerical design choices:
%   * LKF uses Sigma-Point–like propagation, but the points are low-discrepancy
%     Gaussianized lattice nodes rather than deterministic (UKF/CKF) points.
%   * SVD-based square-root extraction is used for numerical stability.
%   * A Cranley–Patterson random shift is applied to the lattice at each step
%     to decorrelate quadrature bias across time.
% ============================================================================

clear; close all; clc;

%% Reproducibility
rng(1,'twister');  % Fix random seed for repeatable results

%% 0) Simulation horizon and dimensions
Tf = 9.0;              % [s] final time
T  = 1e-3;             % [s] sample period
t  = 0:T:Tf;           % time vector
n  = 4;                % # of states
m  = 4;                % # of measurements (here: full-state)

% Pre-allocate data containers
x_true = zeros(n, numel(t));   % true state trajectory
z_meas = zeros(m, numel(t));   % measured outputs
u_sig  = zeros(1, numel(t));   % control/input signal
x_hat  = zeros(n, numel(t));   % LKF state estimate trajectory

%% 1) Noise covariances (tuned to match the stiff nature of state 4)
% Small state noises on x1..x3; large on x4 to emulate hydraulic pressure jitter
Q = diag([1e-12, 1e-10, 1e-9, 1e6]);   % process noise covariance (n x n)
R = 1e3 * Q;                           % measurement noise covariance (m x m)

% Initial covariance for the filter (pessimistic to encourage convergence)
P0 = 10 * Q;

% Generate process and measurement noise sequences
w = mvnrnd(zeros(n,1), Q, numel(t))';   % w_k ~ N(0, Q)
v = mvnrnd(zeros(m,1), R, numel(t))';   % v_k ~ N(0, R)

%% 2) Plant constants (hydraulics – adopted from the original example)
A    = 1.52e-3;    % [m^2]    Piston area
Dp   = 5.57e-7;    % [m^3/rad] Pump displacement
M    = 7.376;      % [kg]     Equivalent mass
V0   = 1.08e-3;    % [m^3]    Hydraulic volume
Beta = 2.07e8;     % [Pa]     Effective bulk modulus
L    = 4.78e-12;   % [m^5/(N*s)] Normal leakage
QL   = 2.41e-6;    % [m^3/s]  Nominal flow rate
a1   = 6.589e4;    % [N*s^2/m^2] Nonlinear friction coeff.
a2   = 2.144e3;    % [N*s/m]     Linear friction coeff.
a3   = 436;        % [N]         Coulomb friction

% Measurement matrix (full-state sensing in this demo)
C = eye(m);

% Input – pump speed profile (square-like), amplitude chosen to excite dynamics.
% We avoid Signal Processing Toolbox and emulate square() by sign(sin(.))
wp = -100 * sign(sin(pi*t));                  % pump speed [rad/s], ±100
u_from = @(x4) (Dp * wp) - sign(x4) * QL;     % leakage-affected input

% Nonlinear plant one-step integrator x_{k+1} = f(x_k, u_k) + w_k
% State ordering: x = [position; velocity; acceleration; pressureDiff]
xNon = @(y) [ ...
    y(1) + T * y(2);  % x1: position
    y(2) + T * y(3);  % x2: velocity
    (1 - T*((a2*V0 + M*Beta*L)/(M*V0))) * y(3) ...
      - T*((A^2*Beta + a2*L*Beta)/(M*V0)) * y(2) ...
      - T/(M*V0) * ( 2*a1*V0*y(2)*y(3) + Beta*L*(a1*y(2)^2 + a3) ) * sign(y(2)) ...
      + T * ((A*Beta)/(M*V0)) * y(5);  % depends on input y(5) = u
    (1/A) * ( a2*y(2) + (a1*y(2)^2 + a3) * sign(y(2)) ) ...
];

%% 3) Simulate the TRUE nonlinear plant (open-loop)
for k = 1:numel(t)-1
    u_sig(k)      = u_from(x_true(4,k));           % leakage-affected input
    x_true(:,k+1) = xNon([x_true(:,k); u_sig(k)]) + w(:,k);  % propagate truth
    z_meas(:,k+1) = C * x_true(:,k+1) + v(:,k+1);            % measure with noise
end

%% 4) Run the Lattice Kalman Filter (LKF)
% Filter configuration
x_hat(:,1) = zeros(n,1);   % initial state estimate
P          = P0;           % initial covariance
N_lattice  = 50;           % # lattice points (typ. 20–200; trade-off cost/accuracy)
a_gen      = 11;           % Korobov generating number (co-prime to N; 7, 8, 9, 11 work)
                           
% Step through time
for k = 1:numel(t)-1
    % One-step LKF update. The LKF internally re-uses xNon, C, Q, R.
    [x_hat(:,k+1), P] = lkf(x_hat(:,k), z_meas(:,k+1), u_sig(k), P, ...
                             xNon, C, Q, R, N_lattice, a_gen);
end

%% 5) Metrics and visualization
err  = x_true - x_hat;
rmse = sqrt(mean(err.^2, 2));

fprintf('\nLattice Kalman Filter RMSE per state:\n');
Ttbl = table((1:n).', rmse, 'VariableNames', {'State','RMSE'});
disp(Ttbl);

% Time series – focus on state 4 (pressure differential), the stiffest state
figure('Color','w'); 
plot(t, x_true(4,:), 'LineWidth', 1.25); hold on;
plot(t, x_hat(4,:),  'k:', 'LineWidth', 1.75);
xlabel('Time [s]'); ylabel('x_4 = \Delta Pressure [Pa]');
title('EHA: True vs. LKF Estimate (state 4)'); legend('True','LKF','Location','best'); grid on;

% Estimation error traces
figure('Color','w');
plot(t, err(4,:), 'LineWidth', 1.25); 
xlabel('Time [s]'); ylabel('Estimation error in x_4 [Pa]');
title('LKF estimation error (state 4)'); grid on;

% Optional: Uncomment to see all states
% figure('Color','w'); tiledlayout(4,1,'TileSpacing','compact');
% for i=1:4
%     nexttile; plot(t,x_true(i,:),'LineWidth',1.0); hold on; plot(t,x_hat(i,:),'k:','LineWidth',1.5);
%     ylabel(sprintf('x_%d',i)); grid on;
% end
% xlabel('Time [s]'); sgtitle('True (solid) vs. LKF (dotted) for all states');

