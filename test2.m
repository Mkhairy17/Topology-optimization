clc
clear all
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
%global stiffness matrix
K=global_matrix2(Lx,Ly,nx,ny);
%%
%Define force vector
F = sparse(2*(ny+1)*(nx+1),1); 
%Define displacement vector
U = sparse(2*(ny+1)*(nx+1),1);
%%
% All degrees of freedom
AllDOF = 1:2*(nx+1)*(ny+1);
%Set Fixed degrees of freedom
FixDOF = [1,2,3,4,5,6];
%Set free degrees of freedom
FreeDOF = setdiff(AllDOF,FixDOF);
%%Define forces
F(9,1)=-10;
%%
%Solve for U
U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
U(FixDOF,:) = 0;
[U_xx,V_yy,Gamaxy]=Calc_str(a,b,nx,ny,U);

%%
%optimiality criteria
ne=nx*ny;
eita1=reshape(full(U_xx),ne,1);
eita2=reshape(full(V_yy),ne,1);
eita12=reshape(full(Gamaxy),ne,1);
P=3;
E=1;
volfrac=0.5;
rhos=0;
l1 = 0; l2 = 1e6;
while (l2-l1) > 1e-9
    lmid = 0.5*(l2+l1);
for ele=1:ne
  rho=(lmid/(P*E*(eita1(ele)^2+4*eita1(ele)*eita12(ele)+4*eita12(ele)^2+2*eita1(ele)*eita2(ele)+4*eita12(ele)*eita2(ele)+eita2(ele)^2)))^(1/(P-1));
  rhom(ele)= rho;
end
rhos=sum(rhom)
if rhos > volfrac*ne
    l2 = lmid;
else
    l1 = lmid; 
end
end



