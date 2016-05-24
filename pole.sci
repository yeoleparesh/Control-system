//Author:- Paresh Yeole
//Dt:-May/2016
//Compute poles of dynamic system
//Syntax : - pole(sys)
//pole(sys) computes the poles p of the SISO or MIMO dynamic system model sys.
//Algorithm:- Used scilab's roots of denominator and converted SS into tf model.
// For MIMO and SISO array cell matrix is created and the poles of the corresponding tf is passed in the cell




function [x]=pole(sys)
     select typeof(sys)
    case "rational" then
        //x=roots(sys.den);
    case "state-space" then
        sys=ss2tf(sys)
       
   
    else
        error(97,1),
    end;
    if(ndims(sys)==3) then
        x=cell(size(sys,'r'),size(sys,'c'),size(sys,3));
    for i=1:size(sys,'r')
     for j=1:size(sys,'c')
         for k=1:size(sys,3)
         x(i,j,k).entries=roots(sys(i,j,k).den);
         end
         
      end
    
    end
    
elseif(ndims(sys)==2) then
    x=cell(size(sys,'r'),size(sys,'c'));
    for i=1:size(sys,'r')
        for j=1:size(sys,'c')
             x(i,j).entries=roots(sys(i,j).den);
            end
        end
      else  
    x=roots(sys.den);
    end
endfunction
