function du=duxyC(U,D,CC)
% compute the emperical derivative of LLC WR to the the elements of Uxy
delta=0.00001;
ll0=LLC(U0,D,CC);

for i=1:length(U0),
  for j=1:length(U0),
    Utmp=U0;
    Utmp(i,j)=Utmp(i,j)+delta;
    llij=LLC(Utmp,D,CC);
    du(i,j)=(llij-ll0)/delta;
  end;
end;
