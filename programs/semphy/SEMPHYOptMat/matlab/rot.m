function r=rot(n,i,j,phi)
% get rotation matrix for one turn

r=eye(n);
r([i,j],[i,j])=[cos(phi),sin(phi);-sin(phi),cos(phi)];
