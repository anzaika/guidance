function r=R(U,D)
% compute R from U,D and pi

PI=U(:,1)'.^2;
pis=diag(sqrt(PI))
pisr=diag(1./(sqrt(PI)))
r=pisr*U'*D*U*pis;