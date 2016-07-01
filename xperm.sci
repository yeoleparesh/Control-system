function [output]=xperm(sys,perm)
    //Function xperm
    //Reorder states in state-space models
    //
    //Calling Sequence:-
    //sys = xperm(sys,Permutation)
    //
    //Parameters:-
    //sys:-
    //      A state-space model of a linear dynamic system.
    //Permuatation:-
    //      The vector Permutation should be a permutation of 1:nx where nx is the number of states in SYS.
    //
    //Description:-
    //  sys = xperm(sys,P) reorders the states of the state-space model sys according to the permutation P.
    // The vector P is a permutation of 1:NX, where NX is the number of states in sys.
    //
    //Examples:-
    //      aa=ssrand(2,3,2);
    //      ss=xperm(aa,[2 1]) 
    //      aa1=ssrand(3,6,8);
    //      ss1=xperm(aa1,[2 5 6 1 4 3 7 8]) 
    //
    //See Also:-
    //       For information about creating state-space models, see ssrand and syslin.
    //
    //Authors
    //  Paresh Yeole emailid:-yeoleparesh@students.vnit.ac.in
    
    
    
    
    
    [lhs,rhs]=argn(0);
    if(rhs<2) then
        error("xperm:Not enough input arguments.");
    end
    select typeof(sys)
    case 'state-space' then
    else
        error("xperm:Undefined function xperm for given type of input arguments ");
    end
    [r,nx]=syssize(sys);
    if ~(or(type(perm)==[1 5 8])) | isequal(gsort(perm),(1:nx))
        error(msprintf(gettext("xperm:The second input argument of the xperm command must be a permutation of (1:%d)."),nx));

    end
    if(length(perm)<>nx) then
        error(msprintf(gettext("xperm:The second input argument of the xperm command must be a permutation of (1:%d)."),nx));
    end

    sys.a = sys.a(perm,perm);
    sys.b = sys.b(perm,:);
    sys.c = sys.c(:,perm);
    output=sys;
endfunction
