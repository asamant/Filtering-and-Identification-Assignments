% Filtering and Identification - final assignment
% Part III - SI State Space Model
%

close all;

load systemMatrices.mat; 
load turbulenceData.mat; % load the data provided

%% Initial

rank_G = rank(G);

% We get the rank of G as 47, whereas the number of columns in G is 49 and
% the number of rows is 72. Hence, G is not full rank.

% Getting the SVD of G:
[Ug, Sg, Vg] = svd(G, 0);

rank_H = rank(H);

% We get the rank of H as 49, and hence we conclude that H is a full rank,
% invertible matrix.

%% Question 5: Sub Identification

% The functions SubId and AOloopSID is present in the same directory as the current
% one.

% Considering the first cell in the given cell array
phiSim_1 = phiIdent{1,19};

% Calculating the dimensions of each phiSim cell
[phi_len, n] = size(phiSim_1);

% Size of G matrix
n_G = size(G,1);

% Number of sample points for phi_sim
T = length(phiSim_1);

% Get sid matrix
sid = G*phiSim_1 + sigmae*randn(n_G,T);

% Choose s and n paramters
s = 25;
n = 23;

[As, Cs, Ks] = SubId(sid, s,n);

%% Question 5: Validation

phiSim_2 = phiSim{1,8};

var_eps = AOloopSID(G,H,As,Cs,Ks,sigmae,phiSim_2);

%% Question 5: Validation over entire dataset

% We are provided with a cell array of 20 datasets for phiSim
num_Datasets = length(phiSim);

variance_closed_loop = 0;
variance_no_control = 0;

for cellIndex = 1:num_Datasets
    phi_currentCell = phiSim{1,cellIndex};    
    
    variance_closed_loop = variance_closed_loop + AOloopRW(G,H, covariance_phi, sigmae, phi_currentCell);
    variance_no_control = variance_no_control + AOloop_nocontrol(phi_currentCell,sigmae,H,G);
end

% Taking the average of the values obtained from all of the provided
% datasets
variance_closed_loop = variance_closed_loop/num_Datasets;
variance_no_control = variance_no_control/num_Datasets;

%% Question 6

% We can clearly see that the variance in the case of the closed-loop
% config is less than a third as the no-control case's, in terms of magnitude

VAF_cumulative = 0;

for cellIndex = 1:num_Datasets
    phi_currentCell = phiSim{1,cellIndex};
    
    VAF_cumulative = VAF_cumulative + VAF_RW(G,H, covariance_phi, sigmae, phi_currentCell);
end

VAF = mean(VAF_cumulative);