function [var_eps] = AOloopRW(G,H, covariance_phi, sigma_eps, phi_sim)
% Variance calculation of an AO system in the closed-loop configuration
% IN
% G     : measurement matrix 
% H     : influence matrix mapping the wavefront on the mirror
% covariance_phi  : covariance matrix of the turbulence wavefront
% sigma_eps : measurement noise parameter for determining its covariance
% phi_sim : simulation data for the wavefront
% OUT
% var_eps : variance of the residual wavefront after taking N_t points
% within the closed-loop operation

% dimension lifted wavefront
n = size(H,1);

% dimension lifted sensor slopes
ns = size(G,1);

% Number of sample points for phi_sim
T = length(phi_sim);

% Residual wavefront
eps_k = zeros(n,T);

u = zeros(n,T);

eps_mean_removed = zeros(n,T); % residual wavefront with mean removed

% Constructing a vector of all the variance values for eps
var_eps = zeros(T,1);

% We have the data for phi_k, and we know the matrix H.
% We also know that eps(k) = phi(k) - Hu(k-1), and our control law requires
% the minimization of eps(k).
% Hence, we follow a linear least-squares approach, minimizing the 2-norm
% of phi(k) - Hu(k-1) to calculate the u(k-1) vector for the sample k.

% An assumption is that H is full rank.



for k = 1:T-1
    % Linear least-squares problem
    u(:,k) = (H'*H)^(-1)*H*phi_sim(:,k+1);

    eps_k(:,k+1) = phi_sim(:,k+1) - H*u(:,k);
    
    eps_mean_removed(:,k+1) = eps_k(:,k+1)-mean(eps_k(:,k+1)); 
%    sk(:,k+1) = G*epsk(:,k+1) + sigmae*randn(ns,1);
    var_eps(k+1) = var(eps_mean_removed(:,k+1));
end

var_eps = mean(var_eps);

end

