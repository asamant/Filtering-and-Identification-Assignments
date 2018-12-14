function [aest,best] = myarx(y,u,n)
%MYARX Takes vectors y and u of the actual outputs and inputs, and an
%"order" n
%   Estimates ARX model coefficients aest and best based on minimizing the
%   difference between predicted and actual output values

[my, ny]=size(y);
[mu,nu]=size(u);
  
if ny>1
    error('y must be a vector.');
end

if nu>1
    error('u must be a vector.') ;
end

if mu ~= my
    error('The vectors u and y must have the same size.');
end

N = my;
N=N(1);
Phi = zeros(N-n, 2*n);

for i=1:N-n
    for j=1:n
        if(i-j > 0)
            Phi(i,j) = y(n+i-j);
            Phi(i, j + n) = u(n+i-j);
        end
    end
end

theta = ((1/N)*(Phi')*Phi)\((1/N)*Phi')*y(n+1:N);

aest = theta(1:n);
best = -theta(n+1:2*n);

end