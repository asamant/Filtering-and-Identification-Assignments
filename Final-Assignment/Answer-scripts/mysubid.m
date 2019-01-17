function [At, Bt, Ct, Dt, x0t, S] = mysubid(y, u, s, n)
%MYSUBID Subspace Identification method

  [N,a1] = size(y);
  [N1,a2] = size(u);
  
  if a1 > 1
    error('The parameter y must have a single column.');
  end
  
  if a2 > 1
    error('The parameter u must have a single column.') ;
  end
  
  if N ~= N1
    error('The vectors u and y must have the same size.');
  end
  
  if n < 1
    error('The order n must be atleast 1.');
  end
  
  if n >= N
    error('The order n must be smaller than the number of measurements.');
  end
  
  if s > N
    error('s must be smaller or equal to the number of measurements.');
  end
  
  if s < n
   error('The order n must be smaller than s.');
  end
  
  Y = zeros(s,N-s+1);
  Uh = zeros(s,N-s+1);
  
  for i = 1:s
    Y(i,:) = y(i:N-s+i)';
    Uh(i,:) = u(i:N-s+i)';
  end
  
  tol = 0.0001;
  if (rank(Uh,tol) < s)
    disp('Warning, inputs do not form a full rank Henkel.');
  end
  
  Proj = eye(N-s+1) - Uh'*inv(Uh*Uh')*Uh;
  
  [U,S,V] = svd(Y*Proj);
  U = U(:,1:n);
  
  Ct = U(1,:);
  At = pinv(U(1:s-1,:))*U(2:s,:);

  [Bt,Dt,x0t] = subidhelp(y,u,At,Ct);
end

