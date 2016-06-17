//author:- Rutuja Moharil & Paresh Yeole  
//function is to find bandwidth of CT system,DT system,SS system for SISO type of systems only.
//[fw]=bandwidth(sys,dbdrop) determines the bandwidth of the system at the given dbdrop. 
//if array of siso system is to be passed then call as:-fw=bandwidth(sys,dbdrop,1)
//"dbdrop must be negative value"
//calling sequence:-
//fw=bandwidth(sys)   ----for SISO transfer function CT,DT and state space --dbdrop=-3db
//fw=bandwidth(sys,dbdrop) --computes bandwidth at given dbdrop
//fw=bandwidth(sys,[],1)   --for SISO array --dbdrop=-3db
//fw=bandwidth(sys,dbdrop,1) ---for SISO array at given dbdrop
function[fw]=bandwidth(sys,varargin)
   //[lhs,rhs]=argn(0);
   n=length(varargin);
   
   if((n>=1 & ( varargin(1)==%nan | varargin(1)>0) ) ) then
       if(n==2 & varargin(1)==[]) then
       else
           error(msprintf(_("dbdrop must be real and negative:bandwidth")));
       end
   
   end
   
   
   if(n==0 | (n==2 & varargin(1)==[]))  then
       
       dbdrop=0.7079;
   //if(varargin(1)==[]) then
       //dbdrop=0.7079;
   elseif(n==1 | n==2) then
      dbdrop=10^(varargin(1)/20);
       end
   select typeof(sys)
case "rational" then
    if(n==2 & varargin(2)==1) then
       o=0;
    else
       o=1;
     end
    
case "state-space" then
        sys=ss2tf(sys)
        o=1;
else
        msprintf("\n")
        error(97,1),
end; 
    
    //bandwidth is defined for only SISO systems.
    if (or(size(sys)<>[1 1]) & o==1) then
      // if ~(typeof(sum(sys)) == "rational") then
        error(msprintf(_("\n %s: Wrong size for input argument #%d: Single input, single output system expected.\n"),"bandwidth",1))
   ///end
   end  
   //[p, m] = size (sys);
     // if (p*m <> 1) then
          //if (typeof(sum(sys)) == "rational") then
       // error(msprintf(_("\n %s: Wrong size for input argument #%d: Single input, single output system expected.\n"),"bandwidth",1))
        //end
   //     end
    
    w=poly(0,'w');
    
    el=1e-7;
    
//niw=((real(horner(sys.num,%i*w)))^2)+((abs(imag(horner(sys.num,%i*w))))^2);
//diw=((real(horner(sys.den,%i*w)))^2)+((abs(imag(horner(sys.den,%i*w))))^2);
if(sys.dt=='c') then
    //t=horner(sys,0);
   if(o==1) then
          try
           t=horner(sys,0);
         catch
            fw=%nan;
             return; 
          end
     
//    if(isinf(t)) then
//    fw=NaN;
//elseif (t==0) then
//    fw= Inf;
//   else 
         niw=(((horner(sys.num,%i*w))))*(conj(horner(sys.num,%i*w)));
         diw=(((horner(sys.den,%i*w))))*(conj(horner(sys.den,%i*w)));
         p=(roots(niw-((t*dbdrop)^2)*diw));
         fw=real(p(find((abs(imag(p))<el)&real(p)>0)));
//end
else
   // aa=size(sys);
    y=ndims(sys);
    //if(y>2)then
 if(y==3) then
    for i=1:size(sys,'r')
        for j=1:size(sys,'c')
          for k=1:size(sys,3)
           try
                t=horner(sys(i,j,k),0);
           catch
                fw(i,j,k)=%nan;
                continue;
             end
 
//                    if(isinf(t)) then
//    fw=NaN;
//elseif (t==0) then
//    fw= Inf;
//    
    //else
            niw=(((horner(sys(i,j,k).num,%i*w))))*(conj(horner(sys(i,j,k).num,%i*w)));
            diw=(((horner(sys(i,j,k).den,%i*w))))*(conj(horner(sys(i,j,k).den,%i*w)));
            p=(roots(niw-((t*dbdrop)^2)*diw));
         k1=real(p(find((abs(imag(p))<el)&real(p)>0)));
    if(k1==[]) then
      fw(i,j,k)=%inf;
   else
    fw(i,j,k)=k1;
   end
//end
end
end
end
//end
else
        for i=1:size(sys,'r')
        for j=1:size(sys,'c')
            
            //for k=1:size(sys,3)
            try
                t=horner(sys(i,j),0);
 
              catch
            
                 fw(i,j)=%nan;
                 continue;
            end

    //if(isinf(t)) then
//    fw=NaN;
//elseif (t==0) then
//    fw= Inf;
//    
//    else
           niw=(((horner(sys(i,j).num,%i*w))))*(conj(horner(sys(i,j).num,%i*w)));
diw=(((horner(sys(i,j).den,%i*w))))*(conj(horner(sys(i,j).den,%i*w)));
p=(roots(niw-((t*dbdrop)^2)*diw));

    k1=real(p(find((abs(imag(p))<el)&real(p)>0)));
    if(k1==[]) then
      fw(i,j)=%inf;
   else
    fw(i,j)=k1;
   end

end
end
end
end
//end
else
    
   if(sys.dt=='d') then
     q=1;
   else
     q=sys.dt
   end 


if(o==1) then
    try

     t=horner(sys,1);
   catch
       fw=%nan;
       end;
    
    //else
  
//z=%e^((%i*w)*q);
//g=horner(sys.num,z)*horner(sys.num,1/z);
//f=horner(sys.den,z)*horner(sys.den,1/z);
//u=roots(((gain*t)^2)*f-g);
//fw=real(u(find((abs(imag(u))<el)&real(u)>0)));
z=poly(0,varn(sys.den));
        sm=simp_mode();
        simp_mode(%f);
        hh=sys*horner(sys,1/z)-(t*dbdrop)^2;
    
        simp_mode(sm);
        //find the numerator roots
        z=roots(hh.num,"e");
        z(abs(abs(z)-1)>el)=[];// retain only roots with modulus equal to 1
        w=log(z)/(%i*q);
        ws=real(w(abs(imag(w))<el&real(w)>0)); //frequency points with unitary modulus
        fw = ws;
else
    for i=1:size(sys,'r')
        for j=1:size(sys,'c')
          for k=1:size(sys,3)
           try
                t=horner(sys(i,j,k),1);
           catch
                fw(i,j,k)=%nan;
                continue;
             end
             z=poly(0,varn(sys(i,j,k).den));
        sm=simp_mode();
        simp_mode(%f);
        hh=sys(i,j,k)*horner(sys(i,j,k),1/z)-(t*dbdrop)^2;
    
        simp_mode(sm);
        //find the numerator roots
        z=roots(hh.num,"e");
        z(abs(abs(z)-1)>el)=[];// retain only roots with modulus equal to 1
        w=log(z)/(%i*q);
        ws=real(w(abs(imag(w))<el&real(w)>0)); //frequency points with unitary modulus
        fw(i,j,k) = ws;
         end
     end
 end
 end
             
        //if fw=[] then it means either the bandwidth is finite or band
//end
end
if(fw==[]) then
    fw=%inf;
    end;
endfunction
    
