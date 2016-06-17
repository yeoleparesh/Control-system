/////Function :-- lyapchol ////////////
////Lyapunov Equation
a  =[-1.6602    1.0973; 1.0973  -2.1947 ];
b=[1.5442 ;0];
r=lyapchol(a,b);
////Generalized Lyapunov Equation
a  =[-1.6602    1.0973; 1.0973  -2.1947 ];
b=[1.5442 ;0];
e=[0.8308    0.5497 ;0.5853    0.9172  ];
r1=lyapchol(a,b,e);
