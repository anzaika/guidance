function R=rotAll(phi)
% get rotation matrix from an angle vector

n=ceil(sqrt(2*length(phi)+1));		% not elegant, but works fine.
R=eye(n);

id=0;
for i=1:(n-1), 
  for j=i+1:n, 
    id=id+1;
    plain(id,:)=[i,j];
  end;
end

for id=1:length(phi),
  r=eye(n);
  r(plain(id,:),plain(id,:))=[cos(phi(id)),sin(phi(id));-sin(phi(id)),cos(phi(id))];
  R=R*r;
end;


