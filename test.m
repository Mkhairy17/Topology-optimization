clc
clear all
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
%global stiffness matrix
K=global_matrix3(Lx,Ly,nx,ny);
%%
%Define force vector
F = sparse(2*(ny+1)*(nx+1),1); 
%Define displacement vector
U = sparse(2*(ny+1)*(nx+1),1);
%%
% All degrees of freedom
AllDOF = 1:2*(nx+1)*(ny+1);
%Set Fixed degrees of freedom
FixDOF = [1,3,5,18];
%Set free degrees of freedom
FreeDOF = setdiff(AllDOF,FixDOF);
%%Define forces
F(2,1)=-1;
%%
%Solve for U
U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
U(FixDOF,:) = 0;
[U_xx,V_yy,Gamaxy]=Calc_str(a,b,nx,ny,U)


