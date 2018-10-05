%% Robust reconstruction of Re=100 cylinder flow
% Comparison of different dictionary/norm combinations
% (not currently used in paper)
%
% Jared Callaham, 2018

load_cylinder

%% Measurement
% Random point measurements
ns = 20;

% Restrict to cylinder wake: final 80% of width, middle 50% of height
sensor_idx = [randperm(round(0.8*ny), ns)' randperm(round(0.5*nx), ns)']; % Choose sensors on restricted area
sensor_idx = [round(0.2*ny)+sensor_idx(:, 1) round(0.25*nx)+sensor_idx(:, 2)];  % Translate to wake
 
% Convert to measurement matrix
C = spdiags(ones(n, 1), 0, n, n);
C = C(sensor_idx(:, 2) + nx*sensor_idx(:, 1), :);  % Sparse measurement matrix

% Plot measurement locations
x = Train(:, 1);
plotCylinder(reshape(x,nx,ny), clim);

%% Alternate dictionary from PCA
[U, Sigma, V] = svd(Train, 'econ');

beta = size(Train,2)/size(Train,1); % Aspect ratio of data matrix
thresh = optimal_SVHT_coef(beta,0) * median(diag(Sigma));

figure(3)
subplot(121)
semilogy(diag(Sigma)); grid on,  hold on
plot([0 size(Sigma, 1)],[thresh thresh],'r--')
subplot(122)
plot(cumsum(diag(Sigma))/sum(diag(Sigma))); grid on, hold on
r = find(diag(Sigma)>thresh, 1, 'last');
plot([r r], [0 1], 'r--')

r = mTrain;
ThetaPCA = U(:, 1:r)*Sigma(1:r, 1:r);

%% K-SVD
% (parameters are hand tuned)
% K = 10;
% nD = 50;
% 
% params.data = Train;
% params.Tdata = K;
% params.dictsize = nD;
% params.iternum = 50;
% params.memusage = 'high';
% [Dksvd,g,~] = ksvd(params,'tr');
% save('ksvd_dict.mat', 'Dksvd', 'nD')

load ksvd_dict.mat

%% Corrupt with noise
eta = 0.25;
x = Test(:, 20);

eps = max(3*eta*rms_vort*sqrt(ns), 0.01); 
noise = eta*rms_vort*randn(size(x));
y = C*(x + noise); % Noisy slice

%% Optimized reconstruction
err = zeros(5, 1);

% Training dictionary, L1 norm
cvx_begin;
    variable s1(mTrain); 
    minimize( norm(s1,1) );
    subject to
        norm(C*Train*s1 - y, 2) <= eps;
cvx_end;
x_hat1 = reconstruct(Train, s1, rms_vort);
err(1) = norm(x_hat1 - x)/norm(x-mean_flow);

% K-SVD dictionary, L1
cvx_begin;
    variable s2(nD);  
    minimize( norm(s2,1) );
    subject to
        norm(C*Dksvd*s2 - y,2) <= eps;
cvx_end;
x_hat2 = reconstruct(Dksvd, s2, rms_vort);
err(2) = norm(x_hat2 - x)/norm(x-mean_flow);


% PCA dictionary, L1
cvx_begin;
    variable s3(r);  
    minimize( norm(s3,1) );
    subject to
        norm(C*ThetaPCA*s3 - y,2) <= eps;
cvx_end;
x_hat3 = reconstruct(ThetaPCA, s3, rms_vort);
err(3) = norm(x_hat3 - x)/norm(x-mean_flow);

% PCA dictionary, L2
cvx_begin;
    variable s4(r);  
    minimize( norm(s4,2) );
    subject to
        norm(C*ThetaPCA*s4 - y,2) <= eps;
cvx_end;
x_hat4 = reconstruct(ThetaPCA, s4, rms_vort);
err(4) = norm(x_hat4 - x)/norm(x-mean_flow);

% Training dictionary, L2
cvx_begin;
    variable s5(mTrain);  
    minimize( norm(s5,2) );
    subject to
        norm(C*Train*s5 - y,2) <= eps;
cvx_end;
x_hat5 = reconstruct(ThetaPCA, s5, rms_vort);
err(5) = norm(x_hat5 - x)/norm(x-mean_flow);

% Print normalized error
err = err/norm(x-mean_flow);
fprintf('Training, L1: \t%0.4f\n', err(1));
fprintf('K-SVD, L1: \t%0.4f\n', err(2));
fprintf('PCA, L1: %0.4f\n', err(3));
fprintf('PCA, L2: \t%0.4f\n', err(4));
fprintf('Training, L2: \t%0.4f\n', err(5));

%% Plots
figure()
%clim = [(min(x)-max(x))/2, (max(x) - min(x))/2];

subplot(231)
plotCylinder(reshape(x+noise,nx,ny), clim); axis off
scatter(sensor_idx(:, 1), sensor_idx(:, 2), 60, 'k', 'filled')
title(sprintf('Test field\nnoise level: %0.2f', eta))

subplot(232)
plotCylinder(reshape(x_hat1, nx, ny), clim); axis off
title(sprintf('Training dictionary, L1\nerr = %0.4f', err(1)));

subplot(233)
plotCylinder(reshape(x_hat2,nx,ny), clim); axis off
title(sprintf('K-SVD, L1\nerr = %0.4f', err(2)))

subplot(234)
plotCylinder(reshape(x_hat3,nx,ny), clim); axis off
title(sprintf('PCA dictionary, L1\nerr = %0.4f', err(3)))

subplot(235)
plotCylinder(reshape(x_hat4, nx, ny), clim); axis off
title(sprintf('PCA dictionary, L2\nerr = %0.4f', err(4)))

subplot(236)
plotCylinder(reshape(x_hat5, nx, ny), clim);
axis off
title(sprintf('Training dictionary, L2\nerr = %0.4f', err(5)))

set(gcf,'Position',[100 100 ny*6 nx*4])
tightfig;
