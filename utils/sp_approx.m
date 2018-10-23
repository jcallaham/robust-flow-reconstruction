function s = sp_approx(y, D, sigma, flow)
%s = SP_APPROX(y, D, eta, avg_energy)
% Sparse approximation y = D*s to some full field x: y=C*x
% INPUTS: D - dictionary (e.g. C*Train)
%         y - measurements
%         sigma - estimate of noise
%         flow - structure containing parameters of flow

 % Estimate of maximum allowable deviation in L2
eps = max(3*sigma*flow.avg_energy*sqrt(length(y))...
    , 0.5*flow.avg_energy*sqrt(length(y)));

m = size(D, 2);

cvx_begin;
    variable s(m);
    minimize( norm(s,1) );
    subject to
        norm(D*s - y, 2) <= eps;
cvx_end;

s = double(s);  % Convert from cvx variable to regular float

end
