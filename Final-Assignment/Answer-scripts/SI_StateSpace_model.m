% Filtering and Identification - final assignment
% Part III - SI State Space Model
%

close all;

load systemMatrices.mat; 
load turbulenceData.mat; % load the data provided

%% Question 5: Sub Identification

% The functions SubId and AOloopSID is present in the same directory as the current
% one.

% Considering the first cell in the given cell array
phiIdent_1 = phiIdent{1,20}; % Trials with different identification datasets

% Calculating the dimensions of each phiSim cell
[phi_len, n] = size(phiIdent_1);

% Size of G matrix
n_G = size(G,1);

% Number of sample points for phi_sim
T = length(phiIdent_1);

% Get sid matrix
sid = G*phiIdent_1 + sigmae*randn(n_G,T);

% Choose s and n paramters
% Chosen using by analysing the SVD using: semilogy(diag(Sigma),'xr');
s = 6;
n = 60;

[As, Cs, Ks] = SubId(sid, s,n); % (Trials with different identification datasets)

%% Question 5: Validation

phiSim_2 = phiSim{1,6}; % Trials with different simulation datasets

[var_eps, VAF] = AOloopSID(G,H,As,Cs,Ks,sigmae,phiSim_2);

%% Question 5: Validation over entire dataset

% We are provided with a cell array of 20 datasets for phiSim
num_Datasets = length(phiSim);

var_eps_total = 0;
VAF_cumulative = 0;

for cellIndex = 1:num_Datasets
    phi_currentCell = phiSim{1,cellIndex};
    [var_currentCell, VAF_currentCell] = AOloopSID(G,H,As,Cs,Ks,sigmae,phi_currentCell);       
    var_eps_total = var_eps_total + var_currentCell;
    VAF_cumulative = VAF_cumulative + VAF_currentCell;
end

% Taking the average of the values obtained from all of the provided
% datasets
fprintf("Results:\n");
Variance_avg = var_eps_total/num_Datasets
VAF_avg = VAF_cumulative/num_Datasets
