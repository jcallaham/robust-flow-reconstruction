function [x_hat, r] = reconstruct(x, Psi, s, flow, rescale)
%[x_hat, r] = RECONSTRUCT(x, Theta, s, flow_params, rescale)
% Calculate normalized residual of the sparse representation of x:
%  x \approx Theta*s
% Note: Uses energy rescaling unless rescale is false

x_hat = Psi*s;  % Sparse reconstruction
if rescale
    x_hat = x_hat*flow.avg_energy/std(x_hat); % Energy rescaling
end
r = norm(x-x_hat)/norm(x+flow.mean_flow); % Normalized L2 residual error

end

