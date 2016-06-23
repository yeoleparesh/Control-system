function[w]=dlyapchol(a,b,e)
//dlyapchol function Compute Cholesky factor of "discrete-time" Lyapunov equations.

//R = dlyapchol(A,B) computes a Cholesky factorization X = R'*R of 
//    the solution X to the Lyapunov matrix equation:
// 
//        A*X*A'- X + B*B' = 0
// 
//    All eigenvalues of A must lie in the open unit disk for R to exist.
// 
//R = dlyapchol(A,B,E) computes a Cholesky factorization X = R'*R of
//    X solving the generalized Lyapunov equation:
// 
//        A*X*A' - E*X*E' + B*B' = 0
// 
//    All generalized eigenvalues of (A,E) must lie in the open unit disk 
//    for R to exist.
//Examples:-
//lyapunov equation
//a  =[-0.25    0.25; 0.6 -0.4];
//b=[1.5442;0];
//r=dlyapchol(a,b)
//generalized lyapunov equation
//a  =[-0.25    0.25; 0.6 -0.4];
//b=[1.5442;0];
//e=[11 22 ;33 44];
//r1=dlyapchol(a,b,e)
//
//

//Author:-Paresh Yeole 
//emailid:-yeoleparesh@students.vnit.in

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
    
    [Gc]=ctr_gram(a,b,'d');
    w=chol(Gc);
case 3 then
    if((size(a)<>size(a')) | (size(e)<>size(a'))) then
        error(msprintf(_("\n A and E matrix must be  square matrices")));
    end
    if((size(a,'r')<>size(b,'r')) | ((size(a,'r')<>size(e,'r')))) then
        error(msprintf(_("\n lyapchol:A , B,E matrices must have same rows")));
    end
    [Gc]=ctr_gram(inv(e)*a,inv(e)*b,'d');
    w=chol(Gc);
else
    error(msprintf(_("\n lyapchol:wrong input arguments")));
end;


endfunction
