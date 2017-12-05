clc
clear all
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
E = 1;
v = 0.3;
%global stiffness matrix
rho_1=ones(ny*nx,1);
rho_2=zeros(ny*nx,1);
P=1;
% All degrees of freedom
AllDOF = 1:2*(nx+1)*(ny+1);
%Set Fixed degrees of freedom
FixDOF = union([1:2:2*(ny+1)],[2*(nx+1)*(ny+1)]);
%Set free degrees of freedom
FreeDOF = setdiff(AllDOF,FixDOF);
%Define force vector
F = sparse(2*(ny+1)*(nx+1),1); 
%%Define forces
F(2,1)=-1;

rho_min = 10^-3;
volfrac=0.5;
Lambda0=300;
%%
while (1)
    K=global_matrix3(Lx,Ly,nx,ny,P,rho_1);
%Define displacement vector
    U = sparse(2*(ny+1)*(nx+1),1);
%Solve for U
    U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
    U(FixDOF,:) = 0;
    strain = Calc_str(a,b,nx,ny,U);
%%
%optimiality criteria
volfractioncalc = @(Lambda) updatedensity(Lambda,rho_1,P,strain,rho_min)/volfrac-1.0;
Lambda=fsolve(volfractioncalc,100);
[volfrac_2,rho_2]=updatedensity2(Lambda,rho_1,P,strain,rho_min);
norm(rho_1-rho_2,'inf')
if norm(rho_1-rho_2,'inf') < 1.0e-1
    break;
end
rho_1=rho_2;
Lambda0=Lambda;
x_new=reshape(rho_2,ny,nx)';
colormap(gray); imagesc(-x_new); axis equal; axis tight; axis off;pause(1e-6);
end


