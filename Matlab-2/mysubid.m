function [At, Bt, Ct, Dt, x0t, S] = mysubid(y, u, s, n)
  
[nmax,a1]=size(y);
[nmax2,a2]=size(u);

%   Initial checks taken from the subidhelp.m file
if a1>1
    error('The parameter y must have a single column.');
end
  
if a2>1
    error('The parameter u must have a single column.') ;
end
  
if nmax2 ~= nmax
    error('The vectors u and y must have the same size.');
end

N=nmax;

y_hankel = zeros(s,N-s+1);
u_hankel = zeros(s,N-s+1);

for i=1:s
    y_hankel(i,:) = y(i:N-s+i)';
    u_hankel(i,:) = u(i:N-s+i)';
end

projection_matrix = eye(N-s+1) - u_hankel'*inv(u_hankel*u_hankel')*u_hankel;

range_matrix = y_hankel*projection_matrix;
S = range_matrix;

[U,S_decomp,V] = svd(y_hankel);

U = U(:,1:n); 

Ct = U(1,:);

At = pinv(U(1:s-1,:))*U(2:s,:);

[Bt,Dt,x0t]=subidhelp(y,u,At,Ct);

end

