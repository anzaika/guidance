function DPHI=dphi(PHI,D,PI,CC)
% compute the emperical derivative WR to the angles

delta=0.00001;
ll0=LL(rotAll(PHI),D,PI,CC);
for i=1:length(PHI),
  p2=PHI;
  p2(i)=p2(i)+delta;
  lli=LL(rotAll(p2),D,PI,CC);
  DPHI(i)=(lli-ll0)/delta;
end;