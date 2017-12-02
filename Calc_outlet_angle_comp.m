function out_angle=Calc_outlet_angle(theta , sigma, a_c ,out_angle_p,indx)
if indx==1
Calc_out = @(out_angle) out_angle-out_angle_p-(0.23*(2*a_c)^2+out_angle/500)*theta/sqrt(sigma);
elseif indx==2
    Calc_out = @(out_angle) out_angle-out_angle_p-(0.23*(2*a_c)^2+out_angle/500)*theta/sigma;
out_angle= fsolve(Calc_out,20);
end