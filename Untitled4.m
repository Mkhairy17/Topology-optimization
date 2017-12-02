clc
clear all
nx=1;
ny=1;
E = 1;
nu = 0.3;
k=[ 1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8); k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3); k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2); k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5);k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4);k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7);k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6);k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];

% All degrees of freedom
AllDOF = 1:2*(nx+1)*(ny+1);
%Set Fixed degrees of freedom
FixDOF = [1 4 7]
%Set free degrees of freedom
FreeDOF = setdiff(AllDOF,FixDOF);
%Define force vector
F = sparse(2*(ny+1)*(nx+1),1); 
%%Define forces
F(8,1)=-1;

U(FreeDOF,:) = KE(FreeDOF,FreeDOF) \ F(FreeDOF,:);
U(FixDOF,:) = 0;