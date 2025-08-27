function [x, P] = lkf(x, z, u, P, xNon, C, Q, R, N, a_gen)
% LKF  — Lattice Kalman Filter (discrete-time, additive Gaussian noise)
% =====================================================================
% [x, P] = lkf(x, z, u, P, xNon, C, Q, R, N, a_gen)
%
% Inputs:
%   x     : (n x 1) prior mean state estimate at time k (x_{k|k-1})
%   z     : (m x 1) measurement at time k
%   u     : scalar (or vector) control input entering xNon via the 5th arg
%   P     : (n x n) prior state covariance at time k (P_{k|k-1})
%   xNon  : function handle, y -> f(y), where y = [x; u] and f returns x_{k+1}
%   C     : (m x n) measurement matrix (z_k = C x_k + v_k)
%   Q     : (n x n) process noise covariance
%   R     : (m x m) measurement noise covariance
%   N     : (1 x 1) number of lattice points (recommend 20–200)
%   a_gen : (1 x 1) Korobov generating number (must be coprime with N)
%
% Outputs:
%   x     : (n x 1) posterior state estimate mean (x_{k|k})
%   P     : (n x n) posterior covariance (P_{k|k})
%
% Algorithm sketch:
%   1. Build Korobov lattice in [0,1)^n: U = (g * j mod N)/N,  j=0..N-1 with g=[1,a,a^2,...].
%   2. Apply a Cranley–Patterson random shift to U to reduce integration bias.
%   3. Map U to Gaussian space via the inverse normal CDF to obtain Y ~ N(0,I).
%   4. Factor P ≈ S S^T (SVD-based square root) and form sigma-like points:
%         X_i = x + S * Y_i,  i=1..N.
%   5. Propagate each X_i through dynamics: X_i <- f([X_i; u]).
%   6. Predicted mean/cov: x^- = mean(X_i),  P^- = Q + (1/N) Σ (X_i - x^-)(... )^T
%   7. Re-linearize sampling around x^- using new P^- and re-sample lattice.
%   8. Form predicted measurements Z_i = C X_i, innovation stats Pzz and Pxz.
%   9. Kalman update: K = Pxz / Pzz,  x = x^- + K(z - z^-),  P = P^- - K Pzz K^T.
%
% Numerical notes:
%   - We symmetrize P before SVD to avoid drift from round-off.
%   - We use an erf-based Gaussian inverse CDF (no toolbox dependency).
%   - We re-sample with a *new* random shift after prediction to reduce bias.
% =====================================================================

% Dimensions
n = size(x,1);
m = size(C,1);

% ----------------------------
% a) Build Gaussian lattice at prior (x, P)
% ----------------------------
% Korobov (rank-1) generator vector g = [1, a, a^2, ... a^{n-1}]^T mod N
g = a_gen .^ (0:n-1)';        % (n x 1)

% U in unit hypercube [0,1)^n, columns are points
jj = 0:(N-1);                 % (1 x N)
U = mod(g * jj, N) / N;       % (n x N)

% Random shift (Cranley–Patterson) to decorrelate bias across time
Shift = rand(1, n);
U = mod1shift(U, Shift);

% Map to Gaussian space via inverse normal CDF (erfinv-based)
Y = sqrt(2) * erfinv(2*U - 1);  % (n x N)  ~ N(0, I)

% Square-root of P by SVD (robust for nearly singular covariances)
P = 0.5*(P + P');     % force symmetry
[U_s, S_s, V_s] = svd(P);
S = 0.5*(U_s + V_s) * sqrt(S_s);    % symmetric "average" square root

% Sigma-like points in Gaussian space
X = x + S * Y;        % (n x N)

% ----------------------------
% b) Prediction through the non-linear dynamics
% ----------------------------
for i = 1:N
    X(:,i) = xNon([X(:,i); u]);     % one-step propagate each sigma-like point
end

x_pred = mean(X, 2);                % predicted mean
P_pred = Q;                         % start with process noise
for i = 1:N
    dx = X(:,i) - x_pred;
    P_pred = P_pred + (1/N) * (dx * dx.');
end

% ----------------------------
% c) Measurement prediction using a *fresh* lattice around (x_pred, P_pred)
% ----------------------------
Shift = rand(1, n);                 % new independent shift
U = mod1shift(U, Shift);            % reuse Korobov structure, new shift
Y = sqrt(2) * erfinv(2*U - 1);

P_pred = 0.5*(P_pred + P_pred.');   % force symmetry
[U_s, S_s, V_s] = svd(P_pred);
S_pred = 0.5*(U_s + V_s) * sqrt(S_s);

X = x_pred + S_pred * Y;            % re-center lattice on predicted mean
Z = C * X;                          % predicted measurement samples
z_pred = mean(Z, 2);

% Innovations and cross-covariances
Pzz = R;
Pxz = zeros(n, m);
for i = 1:N
    dz = Z(:,i) - z_pred;
    dx = X(:,i) - x_pred;
    Pzz = Pzz + (1/N) * (dz * dz.');
    Pxz = Pxz + (1/N) * (dx * dz.');
end

% ----------------------------
% d) Kalman update
% ----------------------------
% Gain (use direct solve; equivalently Pxz * pinv(Pzz) for full-rank Pzz)
K = Pxz / Pzz;

% Posterior
x = x_pred + K * (z - z_pred);
P = P_pred - K * Pzz * K.';
P = 0.5*(P + P.');            % enforce symmetry

end
