function ll=LLC(U,D,CC)
% compute the log likelihood from u,dand data or from R and data
if (nargin==3)
  if (min(size(U))==1)
    A=U;
    U=rotAll(A);
  end  
  Qt = U' * expm(D) * U;
  ll=sum(sum(log(Qt).*CC));
else
  disp('cant do this');
end