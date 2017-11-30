clc
clear all
L=1;
alfa=0.02;
dzeta=0.05;
dt=0.02;
tf=5;
C=alfa*dt*(exp(1)-1)^2/L;
%%
%Initialization of T
num_nodes=L/dzeta+1;
num_times=tf/dt;
zeta=0:dzeta:L;
T=zeros(num_nodes,num_times);
%%
%Initial conditions conditions
T(:,1)=100*sin(0.5*pi.*zeta);
%%
%Boundary conditions
T(1,:)=0;
T(end,:)=100;

%%
for i=2:num_times
    for j=2:num_nodes-1
    T(j,i)=T(j,i-1)+C*exp(-2*zeta(j))*((T(j+1,i-1)-2*T(j,i-1)+T(j-1,i-1))/dzeta^2-(T(j+1,i-1)-T(j-1,i-1)));
    end
end

%Plotting a color map 
pcolor(T)
xlabel('time steps')
ylabel('grid points')