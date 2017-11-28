clc
clear all
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
%global stiffness matrix
rhom=ones(ny,nx)*0.5;
rho_p=zeros(ny,nx);
P=3;

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

E=1;
volfrac=0.5;


while norm(rhom-rho_p,'inf') > 1.0e-3
    
    rho_p=rhom
    K=global_matrix2(Lx,Ly,nx,ny,P,rho_p);
%Define displacement vector
    U = sparse(2*(ny+1)*(nx+1),1);
%Solve for U
    U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
    U(FixDOF,:) = 0;
    strain = Calc_str(a,b,nx,ny,U);
%%
%optimiality criteria

vol = 0;
rho_min=10^-3;

volfractioncalc = @(Lambda) updatedensity(Lambda,strain,P,rho_min,Cm)/volfrac-1.0

fzero(volfractioncalc,Lambda);

l1 = 5; l2 = 1e6;
while (l2-l1) > 1e-101111111
     lmid = 0.5*(l2+l1);

rhos=sum(sum(rhom))
if rhos > volfrac*ne
    l2 = lmid;
else
    l1 = lmid; 
end
end

% rhom_indx1=rhom>1;
% rhom_indx2=rhom<rho_min;
% rhom(rhom_indx1)=1;
% rhom(rhom_indx2)=rho_min;
rhom=reshape(rhom,ny,nx);
end
