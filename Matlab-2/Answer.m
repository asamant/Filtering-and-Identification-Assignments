% SC42025 - FILTERING AND IDENTIFICATION
% Name: ANIKET ASHWIN SAMANT
% Student ID: 4838866

close all;

%% Question 3.1
% Check for full rank of Hankel matrix
n = 6;

N = size(u1);
N = N(1);
u_hankel = zeros(n,N);

for i=1:N-(n-1)
    for j=1:n
        u_hankel(j,i) = u1(i+j-1);
    end
end

rank(u_hankel)

N = size(u2);
N = N(1);
u_hankel = zeros(n,N);

for i=1:N-(n-1)
    for j=1:n
        u_hankel(j,i) = u2(i+j-1);
    end
end

rank(u_hankel)

N = size(u3);
N = N(1);
u_hankel = zeros(n,N);

for i=1:N-(n-1)
    for j=1:n
        u_hankel(j,i) = u3(i+j-1);
    end
end

rank(u_hankel)

N = size(u4);
N = N(1);
u_hankel = zeros(n,N);

for i=1:N-(n-1)
    for j=1:n
        u_hankel(j,i) = u4(i+j-1);
    end
end

rank(u_hankel)
% We see that the rank of the hankel matrix is full(6) only for input
% vectors u1 and u4. Hence, only u1 and u4 are useful for identification.

%% Question 3.2

% s should be greater than or equal to the order of the system so as to
% make the system observable, and should be less than the total number of
% measurements.

%% Question 3.3

% The singular values are used for determining which values of the output
% are not used for system identification.

%% Question 3.4 and 3.5

% We choose s = 7, since it should be greater than the system order, 6.

% Using set 1:

s = 7;
n = 6;

[At1, Bt1, Ct1, Dt1, x0t1, S1] = mysubid(y1, u1, s, n);
[aest1,best1] = myarx(y1,u1,n);

[At4, Bt4, Ct4, Dt4, x0t4, S4] = mysubid(y4, u4, s, n);
[aest4,best4] = myarx(y4,u4,n);

bode(tfse);
hold on;
% a1
sys_id_a1_arx = tf(best1',aest1',-1);
sys_id_a1_subid = ss(At1,Bt1,Ct1,Dt1);
bode(sys_id_a1_arx);
bode(sys_id_a1_subid);

% a4
sys_id_a4_arx = tf(best4',aest4',-1);
sys_id_a4_subid = ss(At4,Bt4,Ct4,Dt4);
bode(sys_id_a4_arx);
bode(sys_id_a4_subid);

% We can see that there is a fair difference between the plots of the ARX
% model and the subspace identification one, the latter being more accurate
% compared with the given transfer function tfse.

%% Question 4

[At, Bt, Ct, Dt, x0t, S] = mysubid(y0, u0, s, 2);
sys_id_a0_subid = ss(At,Bt,Ct,Dt);

% Having an order of 2 seems to give a correct estimate based on the
% transfer function, and hence the order of the system can be said to be 2.
figure(2)
bode(sys_id_a0_subid);