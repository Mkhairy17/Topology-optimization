function [PI_t,tau_t,eta,M1]= Calculate_compressor(R,gama,Cp,m_dot,U1,U2,Tt1,Pt1,A,w_loss_R,w_loss_S,ac_c,ac_s,sigma_R,sigma_s,beta_1p,beta_2p,alfa_1p,alfa_2p,alfa1)
%% Calculate camber angles
theta_c = beta_1p - beta_2p;
theta_s = alfa_1p - alfa_2p;
%% Calculations
% Mass flow parameter equation
Calculate_M = @(M) (m_dot*sqrt(Tt1))/(A*Pt1*cosd(alfa1)) - sqrt(gama/R)*M*(1+M^2*(gama-1)/2)^-((gama+1)/(2*(gama-1)));
M1 = fsolve (Calculate_M , 0.4);
% Calculate  T1 and P1
T1 = Tt1/(1+M1^2*gama*R/(2*Cp));
P1 = Pt1*((1+M1^2*(gama-1)/2))^-(gama/(gama-1));
% Calculate all rotor inlet velocities
C1 = (M1*sqrt(gama*R*T1));
C1x = C1*cosd(alfa1);
C1th = C1*sind(alfa1);
W1th = U1-C1th;
Beta1 = atand(W1th/C1x);
W1 = sqrt(C1x^2+W1th^2);
%% Calculation of relative conditions
M1r = W1/sqrt(gama*R*T1);
Tt1r = T1+W1^2/(2*Cp);
Pt1r = P1/((1+M1r^2*(gama-1)/2))^-(gama/(gama-1));
%% Calculate all rotor outlet velocities % Assumption Cx1=Cx2
Beta2=Calc_outlet_angle(theta_c, sigma_R, ac_c, beta_2p,1);
C2x=C1x;
W2 = C2x/cosd(Beta2);
W2th = C2x*tand(Beta2);
C2th = U2-W2th;
C2 = sqrt(C2th^2+C2x^2);
alfa2 = atand(C2th/C2x);
%% Calculation of rotor outlet temperatures  
Tt2r = Tt1r;
T2 = Tt2r-W2^2/(2*Cp);
Tt2=T2+C2^2/(2*Cp);
%% Mach number calculations at rotor outlet
M2 = C2/sqrt(gama*T2*R);
M2r= W2/sqrt(gama*T2*R);
%% Pressure calculations at the outlet
Pt2r = Pt1r-w_loss_R*(Pt1r-P1);
P2 = Pt2r*((1+M2r^2*(gama-1)/2))^-(gama/(gama-1));
Pt2 = P2/((1+M2^2*(gama-1)/2))^-(gama/(gama-1));
DF_R = 1-(W2/W1)+abs(W1th-W2th)/(2*sigma_R*W1);
%% Stator calculations
% Calculation of rotor outlet velocities
C3x = C2x;
alfa3 = Calc_outlet_angle(theta_s, sigma_s, ac_s, alfa_2p,2);
C3 = C3x/cosd(alfa3);
C3th=C3*sind(alfa3);
%Temperature calculations
Tt3 = Tt2;
T3 = Tt3-C3^2/(2*Cp);
% Mach number calculation
M3 = C3/sqrt(gama*R*T3);
% Pressure calculations
Pt3 = Pt2-w_loss_S*(Pt2-P2);
P3 = Pt3*((1+M3^2*(gama-1)/2))^-(gama/(gama-1));
DF_s = 1-(C3/C2)+abs(C2th-C3th)/(2*sigma_R*C2);
%% Stage performance Calculations
PI_t = Pt3/Pt1;  % Total pressure ratio
PI_s = P3/P1;    % Static pressure ratio
tau_t = Tt3/Tt1;  %Total Temperature ratio
tau_s = T3/T1;    %Static Temperature ratio
eta = ((PI_t^((gama-1)/gama))-1)/(tau_t-1);  % Stage effeciency
end
