function [K]=global_matrix3(Lx,Ly,nx,ny,P,rho)
ne = nx*ny;
row = zeros(64*ne,1);
col = zeros(64*ne,1);
val = zeros(64*ne,1);
indx = 1:64;
%local stiffness matrix 
E = 1;
nu = 0.3;
k=[ 1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8); k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3); k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2); k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5);k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4);k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7);k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6);k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
ie=1;
for j = 1:ny
   for i = 1:nx
       n1 = (ny+1)*(i-1)+j;
       n2 = (ny+1)* i +j;
       edof = [2*n1+1 ; 2*n1+2 ; 2*n2+1 ; 2*n2+2 ; 2*n2-1 ; 2*n2 ; 2*n1-1 ; 2*n1]';
       row(indx) = reshape(repmat(edof,8,1)',64,1);
       col(indx) = reshape(repmat(edof,8,1),64,1);
       val(indx) = (rho(ie)^P)*KE(:);
       indx = indx + 64;
       ie=ie+1;
    end
end    
 
K = sparse(row,col,val);
