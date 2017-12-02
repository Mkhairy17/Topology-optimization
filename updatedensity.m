function [volfrac,rho_2] = updatedensity(Lambda,rho_1,P,strain,rho_min)
E=1;
v=0.3;
Cm = [E*(1-v)/((1+v)*(1-v^2)) v*E*(1-v)/((1+v)*(1-v^2)) 0;v*E*(1-v)/((1+v)*(1-v^2)) E*(1-v)/((1+v)*(1-v^2)) 0;0 0 E/(1+v)];
zeta = 0.1;
eta = 0.5;
ne = size(strain,1);
rho_2 = zeros(ne,1);
for ie = 1:ne
    Strain_energy=strain(ie,:)*Cm*strain(ie,:)';
    b =  P*rho_1(ne)^(P-1)*Strain_energy/Lambda;
    rho_target = b*rho_1(ie);
    if rho_target <= max([(1-zeta)*rho_target rho_min])
        rho_2(ie) = max([(1-zeta)*rho_target rho_min]);
    elseif rho_target >= min([(1+zeta)*rho_target 1])
        rho_2(ie) = min([(1+zeta)*rho_target 1]);
    else
        rho_2(ie) = rho_target;
    end
end

volfrac = sum(rho_2)/ne;

end

