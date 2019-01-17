close all;

load systemMatrices.mat; 
load turbulenceData.mat; % load the data provided

phi_currentCell = phiSim{1,1};
[phi_len, n] = size(phi_currentCell);

% Total number of given datasets with simulated phi values
num_Datasets = length(phiSim);

cov_phi_0_cumulative = 0;
cov_phi_1_cumulative = 0;

% Averaging the covariance values for C_phi(0) and C_phi(1) based on the
% provided datasets
for cellIndex = 1:num_Datasets
    phi_currentCell = phiSim{1,cellIndex};
    
    % Given definition for C_phi(0)
    covariance_phi_0 = zeros(phi_len, phi_len);
    for k = 1:n
        covariance_phi_0 = covariance_phi_0 + (phi_currentCell(:,k)*phi_currentCell(:,k)');
    end

    covariance_phi_0 = covariance_phi_0/n;
    
    cov_phi_0_cumulative = cov_phi_0_cumulative + covariance_phi_0;

    covariance_phi_1 = zeros(phi_len, phi_len);
    
    % Given definition for C_phi(1)
    for k = 2:n
        covariance_phi_1 = covariance_phi_1 + (phi_currentCell(:,k)*phi_currentCell(:,k-1)');
    end

    covariance_phi_1 = covariance_phi_1/(n-1);
    
    cov_phi_1_cumulative = cov_phi_1_cumulative + covariance_phi_1;
    
end

covariance_phi_0 = cov_phi_0_cumulative/num_Datasets;
covariance_phi_1 = cov_phi_1_cumulative/num_Datasets;

%% Question 7

% Computing the values for the Kalman filter
[A, Cw, K] = computeKalmanAR(covariance_phi_0, covariance_phi_1, G, sigmae);

var_KalmanAR = AOloopAR(G,H, covariance_phi_0, sigmae, A, Cw, K, phiSim{1,1});
%% Question 8


