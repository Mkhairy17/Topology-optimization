%Calculate the strains
function [U_xx,V_yy,Gamaxy]=Calc_str(a,b,nx,ny,U)
%derivatives of Interpolation functions with respect to x evaluated at the center of the element 
ye=b/2;
N1_x = ye/(a*b)-(1/a);
N2_x = (1/a)-ye/(a*b);
N3_x = ye/(a*b);
N4_x = -ye/(a*b);
%derivatives of Interpolation functions with respect to y evaluated at the center of the element 
xe=a/2;
N1_y = xe/(a*b)-(1/b);
N2_y = -xe/(a*b);
N3_y = xe/(a*b);
N4_y = (1/b)-xe/(a*b);
% Calculate the strains for each element 
for ely = 1:ny
    for elx = 1:nx
      n1 = (ny+1)*(elx-1)+ely;
      n2 = (ny+1)* elx +ely;
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
      %u_x  is the normal strain in x direction, v_y is the normal strain in y direction 
      u_x = N1_x*Ue(1)+N2_x*Ue(3)+N3_x*Ue(5)+N4_x*Ue(7);
      v_y = N1_y*Ue(2)+N2_y*Ue(4)+N3_y*Ue(6)+N4_y*Ue(8);
      %gama_xy  is the shear strain
      gama_xy = (N1_y*Ue(1)+N2_y*Ue(3)+N3_y*Ue(5)+N4_y*Ue(7))+(N1_x*Ue(2)+N2_x*Ue(4)+N3_x*Ue(6)+N4_x*Ue(8));
      U_xx(ely,elx)=u_x;
      V_yy(ely,elx)=v_y;
      Gamaxy(ely,elx)=gama_xy;
    end
  
end
