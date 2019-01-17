function [A, Cw, K] = computeKalmanAR(C_phi_0, C_phi_1, G, sigma_e)

%COMPUTEKALMANAR Computes parameters for the Kalman filter in the VAR model
%   IN
% C_phi_0 : The auto-covariance matrix for the given phi dataset
% C_phi_1 : The C_phi(1) matrix for the given phi dataset
% G : The G matrix relating epsilon(k) to s(k)
% sigma_e : SD of the noise vector

% OUT
% A : State matrix defining the relation between the current and next
% states
% Cw : Covariance matrix for process noise
% K : Stationary Kalman filter gain 

phi_len = size(C_phi_0, 1); 
A = C_phi_1/C_phi_0;

Cw = (eye(size(A)) - A*A)*C_phi_0;

% For the DARE, we need the Q matrix to be symmetric. So we take the
% average of the Cw matrix and its transpose
Cw = (Cw + Cw')/2;

R = (sigma_e^2)*eye(size(G,1));

E = eye(phi_len);
S = zeros(size(G'));

% Making use of the DARE provided by MATLAB to compute the covariance
% matrix for epsilon
P = dare(A',G',Cw,R,S,E);

% Kalman gain
K = A*P*G'/(G*P*G' + R);

end

