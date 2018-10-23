function [x_corrupt, x_hat, res] = cylinder_reconstruct(x, C, Train, flow, sigma, rho)
%f = CYLINDER_RECONSTRUCT(C, sigma, rho)
%  Solves noisy, corrupt reconstruction problem on the cylinder
%  x - test field
%  C - measurement matrix
%  sigma - noise level
%  rho - % of corrupt measurements
%
% Jared Callaham 2018

%% Measurement over a window in the wake
ns = size(C, 1);
rms_vort = flow.avg_energy;
n = flow.n; mean_flow = flow.mean_flow;

energy_rescale = true;

%% Introduce sparse corruption and dense noise
% Dense, low amplitude noise
noise = sigma*rms_vort*randn(size(x));
x_corrupt = x + noise;

if rho>0
    % Sparse, gross corruption
    corrupt_idx = randperm(n, round(rho*n));
    max_noise = max(max(abs(bsxfun(@plus, Train, mean_flow))));
    corruption = 2*max_noise*(0.5 - rand(size(corrupt_idx)));

    x_corrupt(corrupt_idx) = corruption;
end

% Measure the corrupted data
y = C*x_corrupt;

%% Optimized reconstruction (see Wright et al, 2009)
if rho>0   % Corruption
    B = [C*Train spdiags(ones(ns, 1), 0, ns, ns)];
    w = sp_approx(y, B, sigma, flow);
    [x_hat, res] = reconstruct(x, Train, w(1:flow.mTrain), flow, energy_rescale);
else  % Only dense noise
    s = sp_approx(y, C*Train, sigma, flow);
    [x_hat, res] = reconstruct(x, Train, s, flow, energy_rescale);
end

end