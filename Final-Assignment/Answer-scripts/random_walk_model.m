% Filtering and Identification - final assignment
% Part I - Random Walk Model
%

close all;

load systemMatrices.mat; 
load turbulenceData.mat; % load the data provided

%% Question 1

% We see that we have been provided data corresponding to s0(k) and G, and
% we have to determine the value of phi(k) such that the error, e(k), is
% minimum, given no prior information. This is clearly a linear least
% squares problem.

rank_G = rank(G);

% We get the rank of G as 47, whereas the number of columns in G is 49 and
% the number of rows is 72. Hence, G is not full rank, we conclude.

% Getting the SVD of G:
[U, S, V] = svd(G, 0);

%% Question 5

% In Question 5, we need to know the rank of the system matrix H so as to
% determine a solution for the linear least squares problem minimizing the
% 2-norm of epsilon (k + 1|k)
rank_H = rank(H);

% We get the rank of H as 49, and hence we conclude that H is a full rank,
% invertible matrix.

%% Question 6

% The function AOloopRW is present in the same directory as the current
% one. 

% taking the first cell in phiSim for now
phiSim = phiSim{1,1};

[phi_len, n] = size(phiSim);
covariance_phi = zeros(phi_len, phi_len);

for k = 1:n
    covariance_k = phiSim(:,k)*(phiSim(:,k)');
    covariance_phi = covariance_phi + covariance_k;
end

variance = AOloopRW(G,H, covariance_phi, sigmae, phiSim);