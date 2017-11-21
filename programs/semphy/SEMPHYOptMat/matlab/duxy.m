function du=duxy(U0,D,CC)
% compute the emperical derivative WR to the elements of Ue
delta=0.00001;
ll0=LL(U0,D,CC);

for i=1:length(U0),
  for j=1:length(U0),
    Utmp=U0;
    Utmp(i,j)=Utmp(i,j)+delta;
    llij=LL(Utmp,D,CC);
    du(i,j)=(llij-ll0)/delta;
  end;
end;
