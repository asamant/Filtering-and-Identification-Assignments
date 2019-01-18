% Filtering and Identification - final assignment
% Section 2 - VAR model

close all;

load systemMatrices.mat; 
load turbulenceData.mat; % load the data provided

phi_currentCell = phiSim{1,1};
[phi_len, n] = size(phi_currentCell);

% Total number of given datasets with simulated phi values
num_Datasets = length(phiSim);

%% Question 7

var_KalmanAR = 0;
VAF_Kalman = 0;

for cellIndex = 1:num_Datasets
    phi_currentCell = phiSim{1,cellIndex};
    
    % Given definition for C_phi(0)
    covariance_phi_0 = zeros(phi_len, phi_len);
    for k = 1:n
        covariance_phi_0 = covariance_phi_0 + (phi_currentCell(:,k)*phi_currentCell(:,k)');
    end

    covariance_phi_0 = covariance_phi_0/n;

    covariance_phi_1 = zeros(phi_len, phi_len);
    
    % Given definition for C_phi(1)
    for k = 2:n
        covariance_phi_1 = covariance_phi_1 + (phi_currentCell(:,k)*phi_currentCell(:,k-1)');
    end

    covariance_phi_1 = covariance_phi_1/(n-1);
    
    % Computing the values for the Kalman filter
    [A, Cw, K] = computeKalmanAR(covariance_phi_0, covariance_phi_1, G, sigmae);
    % Computing the variance for the current dataset
    var_KalmanAR = var_KalmanAR + AOloopAR(G,H, covariance_phi_0, sigmae, A, Cw, K, phi_currentCell);
    
    %Computing the VAF values
    VAF_Kalman = VAF_Kalman + VAF_VAR(G,H, sigmae, A, K, phi_currentCell);
end

% Taking the average values over 20 datasets
var_KalmanAR = var_KalmanAR/num_Datasets;
VAF_Kalman = VAF_Kalman/num_Datasets;

%% Question 8

VAF_Kalman
