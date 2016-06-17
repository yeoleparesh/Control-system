function [x]=pole(sys)
    //Author:- Paresh Yeole
//Dt:-May/2016
//Compute poles of dynamic system
//Syntax : - pole(sys)
//pole(sys) computes the poles p of the SISO or MIMO dynamic system model sys.
//For state-space models, the poles are the eigenvalues of the A 
//matrix
//Algorithm:- Used scilab's roots of denominator and converted SS into tf model.
// For MIMO and SISO array cell matrix is created and the poles of the corresponding tf is passed in the cell
//o/p will be a "cell array" 
//in order to access the values convert it into matrix by cell2mat
//Examples:-
//s=poly(0,'s');sysa=syslin('c',(s+2)/(s^2+2*s+3));
//aa=pole(sysa)
//to access the values use:-
//aa=cell2mat(aa);
//aa(1)
//
//sys1=syslin('c',(s+8)/(s^8+3*s+4));
//sysm=[sysa sys1];
//r=pole(sysm);
//r(1)   ---poles of sysm(1,1)
//r(2)   ---poles of sysm(1,2)
//cell2mat(r(2));   ---to disp all the poles 
[lhs,rhs]=argn(0)
    if rhs == 0 then
        error(msprintf(gettext("Not enough input arguments")))
        end

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
//x=cell2mat(x);
endfunction
