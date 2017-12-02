1 %%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE
SIGMUND, OCTOBER 1999 %%%
2 function top(nelx,nely,volfrac,penal,rmin);
3 % INITIALIZE
4 x(1:nely,1:nelx) = volfrac;
5 loop = 0;
6 change = 1.;
7 % START ITERATION
8 while change > 0.01
9 loop = loop + 1;
10 xold = x;
11 % FE-ANALYSIS
12 [U]=FE(nelx,nely,x,penal);
13 % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
14 [KE] = lk;
15 c = 0.;
16 for ely = 1:nely
17 for elx = 1:nelx
18 n1 = (nely+1)*(elx-1)+ely;
19 n2 = (nely+1)* elx +ely;
20 Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;
2*n2+2; 2*n1+1;2*n1+2],1);
21 c = c + x(ely,elx)^penal*Ue’*KE*Ue;
22 dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*
Ue’*KE*Ue;
23 end
24 end
25 % FILTERING OF SENSITIVITIES
26 [dc] = check(nelx,nely,rmin,x,dc);
27 % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
28 [x] = OC(nelx,nely,x,volfrac,dc);
29 % PRINT RESULTS
30 change = max(max(abs(x-xold)));
31 disp([’ It.: ’ sprintf(’%4i’,loop) ’ Obj.: ’
sprintf(’%10.4f’,c) ...
32 ’ Vol.: ’ sprintf(’%6.3f’,sum(sum(x))/
(nelx*nely)) ...
33 ’ ch.: ’ sprintf(’%6.3f’,change )])
34 % PLOT DENSITIES
35 colormap(gray); imagesc(-x); axis equal; axis
tight; axis off;pause(1e-6);
36 end
37 %%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%
38 function [xnew]=OC(nelx,nely,x,volfrac,dc)
39 l1 = 0; l2 = 100000; move = 0.2;
40 while (l2-l1 > 1e-4)
41 lmid = 0.5*(l2+l1);
42 xnew = max(0.001,max(x-move,min(1.,min(x+move,x.
*sqrt(-dc./lmid)))));
43 if sum(sum(xnew)) - volfrac*nelx*nely > 0;
44 l1 = lmid;
45 else
46 l2 = lmid;
47 end
48 end
49 %%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%
50 function [dcn]=check(nelx,nely,rmin,x,dc)
51 dcn=zeros(nely,nelx);
52 for i = 1:nelx
53 for j = 1:nely
54 sum=0.0;
55 for k = max(i-round(rmin),1):
min(i+round(rmin),nelx)
56 for l = max(j-round(rmin),1):
min(j+round(rmin), nely)
57 fac = rmin-sqrt((i-k)^2+(j-l)^2);
58 sum = sum+max(0,fac);
59 dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)
*dc(l,k);
60 end
61 end
62 dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
63 end
64 end
65 %%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%
66 function [U]=FE(nelx,nely,x,penal)
67 [KE] = lk;
68 K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*
(nely+1));
69 F = sparse(2*(nely+1)*(nelx+1),1); U =
sparse(2*(nely+1)*(nelx+1),1);
70 for ely = 1:nely
71 for elx = 1:nelx
72 n1 = (nely+1)*(elx-1)+ely;
73 n2 = (nely+1)* elx +ely;
74 edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1;
2*n2+2;2*n1+1; 2*n1+2];
75 K(edof,edof) = K(edof,edof) +
x(ely,elx)^penal*KE;
76 end
77 end
78 % DEFINE LOADSAND SUPPORTS(HALF MBB-BEAM)
79 F(2,1) = -1;
80 fixeddofs = union([1:2:2*(nely+1)],
[2*(nelx+1)*(nely+1)]);
81 alldofs = [1:2*(nely+1)*(nelx+1)];
82 freedofs = setdiff(alldofs,fixeddofs);
83 % SOLVING
127
84 U(freedofs,:) = K(freedofs,freedofs) \
F(freedofs,:);
85 U(fixeddofs,:)= 0;
86 %%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%
87 function [KE]=lk
88 E = 1.;
89 nu = 0.3;
90 k=[ 1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...
91 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
92 KE = E/(1-nu^2)*
[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
93 k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
94 k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
95 k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
96 k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
97 k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
98 k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
99 k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];