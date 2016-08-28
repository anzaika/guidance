function ll=LL(U,D,T,CC)
% compute the log likelihood from u,d,t and data or from R and data
if (nargin==4)
  if (min(size(U))==1)
    A=U;
    U=rotAll(A);
  end  
  PI=U(1,:).^2;
  ll=sum(sum(log(diag(PI)*expm(R(U,D)*T)).*CC));%/sum(sum(CC));
else
  CC=D;
  r=U;
  PI=mean(expm(r*1000));
  ll=sum(sum(log(diag(PI)*expm(r*T)).*CC));%/sum(sum(CC));
end