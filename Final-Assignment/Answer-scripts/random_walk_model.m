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
% the number of rows is 72. Hence, G is clearly not full rank.


