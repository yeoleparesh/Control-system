function[w]=lyapchol(a,b,e)
//Calling Sequence
//  R=lyapchol(A,B)
//  R=lyapchol(A,B,E)
//Description
//lyapchol function Compute Cholesky factor of "continuous-time" Lyapunov equations.
//
//R=lyapchol(A,B)
//A X  +  XA'  +  B B'  =  0           (Lyapunov Equation)
//X=R*R'
//
//R=lyapchol(A,B,E)
//A X E'  +  E XA'  +  B B'  =  0      (Generalized Lyapunov Equation)
//X=R*R'
//Examples
//
// 1.a  =[- 1.6602    1.0973; 1.0973  - 2.1947 ];
//  b=[1.5442 ;0];
//  r=lyapchol(a,b);
//  r=[0.9591402    0.3183234 ;0.           0.2265428  ];
//
//
// 2.a  =[- 1.6602    1.0973; 1.0973  - 2.1947 ];
//   b=[1.5442 ;0];
//   e=[0.8308    0.5497 ;0.5853    0.9172  ];
//   r=lyapchol(a,b,e);
//   r=[1.00442  - 0.0600239  ;0.         0.5714961  ];
//See also
// dlyap,dlyapchol
//Authors
 //Paresh Yeole
[lhs rhs]=argn(0);
    if rhs<2  then
        error(msprintf(gettext("Not enough input arguments")));
    end
    
select argn(2)
case 2 then
    if(size(a)<>size(a')) then
        error(msprintf(_("\n A matrix must be a square matrix")));
    end
    if(size(a,'r')<>size(b,'r')) then
        error(msprintf(_("\n lyapchol:A and B matrices must have same rows")));
    end
    
    [Gc]=ctr_gram(a,b);
    w=chol(Gc);
case 3 then
    if((size(a)<>size(a')) | (size(e)<>size(a'))) then
        error(msprintf(_("\n A and E matrix must be  square matrices")));
    end
    if((size(a,'r')<>size(b,'r')) | ((size(a,'r')<>size(e,'r')))) then
        error(msprintf(_("\n lyapchol:A , B,E matrices must have same rows")));
    end
    [Gc]=ctr_gram(inv(e)*a,inv(e)*b);
    w=chol(Gc);
else
    error(msprintf(_("\n lyapchol:wrong input arguments")));
end;


endfunction



 
