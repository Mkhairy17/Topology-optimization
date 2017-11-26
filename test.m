clc
clear all
Lx=input('length of x=');
Ly=input('length of y=');
nx=input('nx=');
ny=input('ny=');
%global stiffness matrix
K=global_matrix2(Lx,Ly,nx,ny);
%convert from sparse matrix to full matrix
% FullM=full(K);

