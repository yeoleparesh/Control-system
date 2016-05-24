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
        if(ndims(p)==3) then
            x=cell(size(p,'r'),size(p,'c'),size(p,3));
            for i=1:size(p,'r')
             for j=1:size(p,'c')
              for k=1:size(p,3)
                  try
                  x(i,j,k).entries=horner(p(i,j,k),0);
              catch
                  x(i,j,k)=%inf;
                  end
                  end
             end
             
            end
            else
            try
    x=horner(p,0)
catch
    x=%inf;
end
end
else
    if(ndims(p)==3) then
            x=cell(size(p,'r'),size(p,'c'),size(p,3));
            for i=1:size(p,'r')
             for j=1:size(p,'c')
              for k=1:size(p,3)
                  try
                  x(i,j,k).entries=horner(p(i,j,k),1);
              catch
                  x(i,j,k)=%inf;
                  end
                  end
             end
             
            end
    
        else
            //if system is discrete then put z=1 for dc gain
  try
    x=horner(p,1)
catch
    x=%inf;
    end
    end
end

endfunction
