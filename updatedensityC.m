function [volfrac,rho_2] = updatedensityC(Lambda,rho_1,P,strain,rho_min)
E=1;
v=0.3;
Cm = (E/1-v^2)*[1 v 0;v 1 0;0 0 1-v];
zeta = 0.05;
eta = 1;
ne = size(strain,1);
rho_2 = zeros(ne,1);
for ie = 1:ne
    Strain_energy=rho_1(ie)^(2*P)*(strain(ie,:)*Cm*strain(ie,:)');
    rho_target = (P*Strain_energy/Lambda)^(1/(P+1));
    rho_2(ie)= min([max([rho_target rho_min]) 1]);
end

volfrac = sum(rho_2)/ne;

end

