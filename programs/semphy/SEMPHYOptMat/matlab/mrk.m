function R=mrk(phi,k,v)
% compute M^k_v for v='cos'|'sin'|'one'

n=ceil(sqrt(2*length(phi)+1));
R=eye(n);

id=0;
for i=1:(n-1), 
  for j=i+1:n, 
    id=id+1;
    plain(id,:)=[i,j];
  end;
end

if (v=='cos')  
  drot=zeros(n,n);
  drot(plain(k,:),plain(k,:))=[1,0;0,1];
else if (v=='sin')
  drot=zeros(n,n);
  drot(plain(k,:),plain(k,:))=[0,1;-1,0];
else if (v=='one')
  drot=eye(n);
  drot(plain(k,:),plain(k,:))=[0,0;0,0];
else 
  disp('v is unknown, exitig');
  return
end;
end;
end;


for id=1:length(phi),
  r=eye(n);
  if (id == k)
    r=drot;
  else
    r(plain(id,:),plain(id,:))=[cos(phi(id)),sin(phi(id));-sin(phi(id)),cos(phi(id))];
  end;
  R=r*R;
end;


