function [volfrac,rho_p] = updatedensity(Lambda,strain,P,rho_min,Cm) 
zeta = 0.1;
eta = 0.5;
ne = size(strain,1);
rho_p = zeros(ne,1);
for ie = 1:ne
    b =  P*E*(rhom(ie)^(P-1))*(eita1(ie)^2+4*eita1(ele)*eita12(ele)+4*eita12(ele)^2+2*eita1(ele)*eita2(ele)+4*eita12(ele)*eita2(ele)+eita2(ele)^2)/Lambda;
    rho_target = b^eta*rhom(ie);
    rho_u = min([(1+zeta)*rho_target 1]);
    rho_l = max([(1-zeta)*rho_target rho_min]);    
    rho_p(ie) = max(min(rho_target,rho_u),rho_l);
end

volfrac = sum(rho_p)/ne;

end

