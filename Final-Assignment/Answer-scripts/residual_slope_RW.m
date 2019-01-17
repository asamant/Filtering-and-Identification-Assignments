% Part 4 - Critical thinking

close all;

load systemMatrices.mat; 
load turbulenceData.mat; % load the data provided

%[U, S, V] = svd(G*H, 0);

[U_G, S_G, V_G] = svd(G);

[U, S, V] = svd(G*H);
U1 = U(:,1:47);
V1 = V(:,1:47);
Sigma_SVD = S(1:47,1:47);

% dimension lifted wavefront
n_H = size(H,1);

% dimension lifted sensor slopes
n_G = size(G,1);

phi_sim = phiSim{1,1};

phi_size = size(phi_sim,1);

covariance_phi = zeros(phi_size, phi_size);

% Number of sample points for phi_sim
T = length(phi_sim);

u = zeros(n_H,T);
for k = 1:T
    covariance_phi = covariance_phi + (phi_sim(:,k)*phi_sim(:,k)');
end

covariance_phi = covariance_phi/T;

% This term is multiplied with the slope vector to obtain the eps_hat(k|k)
% value

s_k = zeros(n_G,length(phi_sim));
eps_k = zeros(n_H,length(phi_sim));

s_k(:,1) = G*phi_sim(:,1) + sigmae*randn(n_G,1);

eps_pred_multiplier_matrix = (covariance_phi*G'/(G*covariance_phi*G' + (sigmae^2)*eye(n_G)));

% Linear least-squares solution 
u(:,1) = (V1/Sigma_SVD)*U1'*(G*eps_pred_multiplier_matrix*s_k(:,1));

var_s = zeros(T,1);

for k = 2:T
    s_k(:,k) = G*(phi_sim(:,k) - H*u(:,k-1)) + sigmae*randn(n_G,1);
    Y_k = G*eps_pred_multiplier_matrix*s_k(:,k) + G*H*u(:,k-1);
    u(:,k) = (V1/Sigma_SVD)*U1'*Y_k;
    eps_k(:,k) = phi_sim(:,k) - H*u(:,k-1);
    
    var_s(k) = var(eps_k(:,k)) - mean(eps_k(:,k));
end

var_s = mean(var_s);
