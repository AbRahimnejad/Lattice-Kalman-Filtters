% =====================================================================
% Helper: Cranley–Patterson shift modulo 1 (vectorized)
% =====================================================================
function Ps = mod1shift(P, delta)
% P: (n x N) lattice points in [0,1)
% delta: (1 x n) shift vector (row)
delta = reshape(delta, size(P,1), 1);      % make it (n x 1)
delta = repmat(delta, 1, size(P,2));       % broadcast to (n x N)
Ps = mod(P + delta, 1);                    % wrap into [0,1)
end
