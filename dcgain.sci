function[x]=dcgain(p)
    //function is to find the dc gain of LTI system
//k = dcgain(sys) computes the DC gain k of the LTI model sys.
//DC gain is infinite for systems with integrators.
//author:- Paresh Yeole emailid:-yeoleparesh@students.vnit.ac.in
[lhs rhs]=argn(0);
if rhs<1 then
    error("dcgain:input parameter as a dynamic system is expected");
end
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
            //x=cell(size(p,'r'),size(p,'c'),size(p,3));
            //t=cell(size(p,'r'),size(p,'c'),size(p,3));
            for i=1:size(p,'r')
             for j=1:size(p,'c')
              for k=1:size(p,3)
                  try
                   t=horner(p(i,j,k),0);
              catch
                   //x(i,j,k).entries=1;
                   //continue;
                   
                   x(i,j,k)=%inf;
                   //disp(t);
                  continue;
 //return;
                  end
          
          x(i,j,k)=t;
          
             //=y;
//             if(t==1) then
//                 x(i,j,k).entries=%inf;
//                 end
              
                  end
             end
             
            end
            
       elseif(ndims(p)==2) then
            //x=cell(size(p,'r'),size(p,'c'));
            for i=1:size(p,'r')
            for j=1:size(p,'c')
             
                 try
                  t=horner(p(i,j),0);
             catch
                  x=%inf;
                  continue;
                end
                  x(i,j)=t
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
    //discrete systems
    if(ndims(p)==3) then
           // x=cell(size(p,'r'),size(p,'c'),size(p,3));
            for i=1:size(p,'r')
             for j=1:size(p,'c')
              for k=1:size(p,3)
                  try
                  t=horner(p(i,j,k),1);
              catch
                  x(i,j,k)=%inf;
                  continue;
              end
              x(i,j,k)=t;
                  end
             end
             
            end
   elseif(ndims(p)==2) then
        //x=cell(size(p,'r'),size(p,'c'));
            for i=1:size(p,'r')
             for j=1:size(p,'c')
             
                  try
                  t=horner(p(i,j),1);
              catch
                  x(i,j)=%inf;
                  continue;
              end
              x(i,j,k)=t;
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
