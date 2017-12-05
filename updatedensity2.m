function [volfrac,rho_2] = updatedensity2(Lambda,rho_1,P,strain,rho_min)
E=1;
v=0.3;
Cm = (E/1-v^2)*[1 v 0;v 1 0;0 0 1-v];
ne = size(strain,1);
rho_2 = zeros(ne,1);
for ie = 1:ne
    Strain_energy=strain(ie,:)*Cm*strain(ie,:)';
    rho_2(ie)= (Lambda/(Strain_energy*P))^(1/(P-1));
end

volfrac = sum(rho_2)/ne;

end

