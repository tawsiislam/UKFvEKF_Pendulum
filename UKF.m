function [mu,Sigma] = UKF(mu_prev,Sigma_prev,z,R,Q,m,L,dt)
ndim = length(mu_prev(:,1));   % Number of states we are estimating,, giving the dimension
alpha = 1;    % Changable param, decides distance from mean along eigenvector, must be [0,1)
kappa = 0;    % Changable param, same property as alpha, must be non-negative
beta = 2; % Default value, should always be beta = 2

lamda = alpha^2*(ndim+kappa)-ndim;

% Create sigma points
% X = zeros(ndim,5); % Contains sigma points
gamma = sqrt(ndim+lamda);
Sigma_chol = gamma*chol(Sigma_prev,"lower"); % Square root of cov_mat is found by Cholesky factorization
X = [mu_prev mu_prev+Sigma_chol mu_prev-Sigma_chol];

% Setting the weights
Wm = zeros(1,2*ndim+1); % Contains weights used for calculating mean
Wc = zeros(1,2*ndim+1);    % Contains weights used for calculating covariance matrix
Wm(1) = lamda/(ndim+lamda);
Wc(1) = lamda/(ndim+lamda) + (1-alpha^2+beta);
Wm(2:end) = 1/(2*(ndim+lamda));
Wc(2:end) = 1/(2*(ndim+lamda));

W_cov = Wc;

% Get "predicted state"
X_barprim=state_trans(X, m, L, dt) + 1*normrnd(0,0.1,size(mu_prev,1),1); % propagate Sigma points through state_trans
mu_bar= X_barprim*Wm';

% Estimated measurement
if length(L) == 1
    C = [1 0];
elseif length(L) == 2
    C = [1 0 0 0;
         0 1 0 0];
end
Z_bar = C*X_barprim;
Z_hat = Z_bar*Wm';

P_ = R; Pyy = Q; Pxy = zeros(ndim,1);  
for i_iters = 1:2*ndim+1 
    M = X_barprim(:,i_iters) - mu_bar;
    N = Z_bar(:,i_iters) - Z_hat;
   
   P_ = P_ + W_cov(i_iters)*(M*M'); % Predicted covariance
   Pyy = Pyy + W_cov(i_iters)*(N*N'); % Measurement/innovation covariance
   Pxy = Pxy + W_cov(i_iters)*(M*N'); % State-measurement cross-covariance
end

Kalman = Pxy/Pyy;
innov = z-Z_hat;
mu = mu_bar+Kalman*innov;
Sigma = P_ - Kalman*Pyy*Kalman';

end