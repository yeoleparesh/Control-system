function X = dlyap(A,B,C,E)
//    DLYAP  Solve discrete Lyapunov equations.
//
//   X = DLYAP(A,C) solves the discrete Lyapunov matrix equation:
//
//       A'*X*A-X +C=0
//
//   X = DLYAP(A,B,C) solves the Sylvester equation:
//
//       A*X*B - X +C =0 
//
//   X = DLYAP(A,C,[],E) solves the generalized discrete Lyapunov equation:
//
//       A*X*A' - E*X*E'  +C =0 
//
//Note:-The discrete-time Lyapunov equation has a unique solution X, for any C = C'
//      if and only if λi(A)λj(A)<>1, for i,j = 1, . . . , n.
//
//References:-
//          http://stanford.edu/class/ee363/notes/lq-lyap-notes.pdf
//
////Examples:-
//A = [3.0   1.0   1.0
//      1.0   3.0   0.0
//      0.0   0.0   3.0];
//
// B = [25.0  24.0  15.0
//      24.0  32.0   8.0
//      15.0   8.0  40.0];
//
// X = dlyap (A, B);
//
//
// Sylvester
// A = [1.0   2.0   3.0 
//      6.0   7.0   8.0 
//      9.0   2.0   3.0];
//
// B = [7.0   2.0   3.0 
//      2.0   1.0   2.0 
//      3.0   4.0   1.0];
//
// C = [271.0   135.0   147.0
//      923.0   494.0   482.0
//      578.0   383.0   287.0];
//
// X = dlyap (A, B, C);
//
////Author:Paresh Yeole emailid:-yeoleparesh@students.vnit.ac.in
    [lhs rhs]=argn(0);
    select rhs
    case 2 then
        X=lyap(A.',-B,'d');
    case 3 then
        X=sylv(-A,B,C,'d');
    case 4 then
        if C<>[] then
            error("dlyap:for generalised lyapunov equation the third input argument must be empty square brackets []");
        end
        
        X=lyap(inv(E)*(A.'),inv(E)*(-B),'d');
    else
        error("dlyap:wrong no. of input arguments");
    end
    
       //X=linmeq(2,A,C,flag)
endfunction
