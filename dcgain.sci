//author:- Paresh Yeole
//function is to find the dc gain of LTI system
//k = dcgain(sys) computes the DC gain k of the LTI model sys.
//DC gain is infinite for systems with integrators.
function[x]=dcgain(p)
    select typeof(p)
    case "rational" then
    case "state-space" then
        p=ss2tf(p)
    else
        msprintf("\n")
        error(97,1),
    end;
    if(p.dt=='c') then
        //if system is continuous then put s=o for dc gain
    x=horner(p,0)
else
    //if system is discrete then put z=1 for dc gain
    x=horner(p,1)
end

endfunction
