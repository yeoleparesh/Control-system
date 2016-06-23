// Function dlyapchol
//lyapunov equation
a  =[-0.25    0.25; 0.6 -0.4];
b=[1.5442;0];
r=dlyapchol(a,b)
//generalized lyapunov equation
a  =[-0.25    0.25; 0.6 -0.4];
b=[1.5442;0];
e=[11 22 ;33 44];
r1=dlyapchol(a,b,e)
