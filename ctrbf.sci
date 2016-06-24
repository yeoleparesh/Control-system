function [aout,bout,cout,t,k] = ctrbf(a, b, c, tol)
    //function ctrbf returns the Controllability staircase form of the system having a,b,c matrices
    //
    //The approach to find whether the system is controllable from its cotrollability matrix is not reliable
    //especially if the eigenvalues are sensitive to small perturbations. Therefore , for more robust approach
    //we use Controllability staircase form i.e. the staircase reduction of (A,B) pair.
    // Here rank is revealed directly from the submatrices in the staircase form.
    //
    //
    //Calling sequence:-
    //[aout,bout,cout,t,k] = ctrbf(a, b, c)
    //                      returns the decomposition of matrices a,b,c into 
    //                      controllable/uncontrollable subspaces
    //
    //[aout,bout,cout,t,k] = ctrbf(a, b, c, tol)
    //                       Uses the tolerance tol for returning the decomposition
    //
    //
    //The last output K is a vector of length n containing the number of controllable states 
    //The number of controllable states is SUM(K).
    ///
    //
    //Reference:-
    //  https://www8.cs.umu.se/~stefanj/pub/SJohansson05.pdf
    //
    //
    //Examples:-
    //A=[1 1;4 -2];B=[1 -1;1 -1];C=[1 0;0 1];
    //[a b c k]=ctrbf(A,B,C)
    //
    //[a b c k t]=ctrbf(A,B,C,10)
    //
    //Author:-Paresh Yeole 
    //emailid:-yeoleparesh@students.vnit.ac.in
    
    
    
    
    
    
       
    [lhs rhs]=argn(0);
    [rowsa,colsa] = size(a);
    [rowsb,colsb] = size(b);
    mat = eye(rowsa,rowsa);
    ajn1 = a;
    bjn1 = b;
    rojn1 = colsb;
    del1 = 0;
    sig1 = rowsa;
    k = zeros(1,rowsa);
    if rhs == 3
        tol = rowsa*norm(a,1)*%eps;
    end
    for jj = 1 : rowsa
        [uj,sj,vj] = svd(bjn1);
        [rsj,csj] = size(sj);
        p = flipdim((eye(rsj,rsj))',1);
        uj = uj*p;
        bb = uj'*bjn1;
        roj = rank(bb,tol);
        [rbb,cbb] = size(bb);
        Rrank = rbb - roj;
        sig1 = Rrank;
        k(jj) = roj;
        if roj == 0 
            break; 
        end
        if Rrank == 0
            break;
        end
        abxy = uj' * ajn1 * uj;
        aj   = abxy(1:Rrank,1:Rrank);
        bj   = abxy(1:Rrank,Rrank+1:Rrank+roj);
        ajn1 = aj;
        bjn1 = bj;
        [ruj,cuj] = size(uj);
        ptj = mat *[uj zeros(ruj,del1);zeros(del1,cuj) eye(del1,del1)];
        mat = ptj;
        deltaj = del1 + roj;
        del1 = deltaj;
    end

    t = mat';
    aout = t * a * t';
    bout = t * b;
    cout = c * t';
endfunction
