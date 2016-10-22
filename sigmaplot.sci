function H=sigmaplot(varargin)
//Calling Sequence
// sigmaplot(sys)
// sigmaplot(sys,{Wmin,Wmax})
// sigmaplot(sys,W)
// sigmaplot(sys1,sys2,sys3,sys4,...,W)
//
//Parameters
// sys - dynamic system
// Wmin - Minimum frequency in radians/unit from which the singular values of the system are to be calculated.
// Wmax - Maximum frequency in radians/unit upto which the singular values of the system are to be calculated.
// W - Range of frequencies in which the singular values of the system has to be calculated.
// V - Singular values of the systems in the range of frequency W.
//
//Description
// sigmaplot function plots or calculates Singular values of dynamic systems.
// When frequency range is not specified, it is chosen automatically.
//
//Examples
// s=%s;
// sys=syslin('c',1/(1+s));
// sigmaplot(sys)
// sigmaplot(sys,sys^2)
// sigmaplot(sys,20:40)
//
//See also
// sigma,bode.
//
//Authors
// Paresh Yeole
// email-id:- yeoleparesh@students.vnit.ac.in

 [lhs,rhs]=argn(0);
     x=logspace(-3,3);
    //________________if freq vector or freqfinal is given//__________________//________________//
    if rhs > 1 then
        if typeof(varargin($)) == 'constant' then
            if ~isreal(varargin($)) then
                error("sigma:The freq values must be real and positive")
            end
            if size(varargin($)) == [1 2] then
                if varargin($) <= 0 then
                    error(msprintf(gettext("%s: The final time value must be a positive real number.\n"),"sigma"))
                end
                if(varargin($)(1)>=varargin($)(2)) then
                    error("sigma:The frequency interval must be specified as {WMIN,WMAX} where WMIN and WMAX are positive frequencies satisfying WMIN<WMAX.");
                end
                
                x=varargin($)(1):0.01:varargin($)(2);
           elseif (isequal(size(varargin($)),[1 2]) == %f & size(varargin($),'*')<>1) then
               x=varargin($);
           elseif (size(varargin($),'*')==1) then
                error("sigma:The frequency interval must be specified as {WMIN,WMAX} where WMIN and WMAX are positive frequencies satisfying WMIN<WMAX.");
                    
              end
            end
        end
    //______________________________________________________________//
    
    
    
    
    CIindex=1;
   
//____________PLOTS__________________________________//
    if lhs==1 then
        for i=1:rhs
            CIindex=CIindex+1;
            if or(typeof(varargin(i))==['string','boolean']) then
                continue;
            end
            
            //__________________SiSo system & MIMO systems__________________________//
            if(or(typeof(varargin(i))==['rational','state-space']) & (i==rhs | typeof(varargin(i+1))<>'boolean' )) then
                y=svplot(varargin(i),x);
                y=20*log10(y);
                if (i<>rhs & typeof(varargin(i+1))=='string') then
                    plot(x,y,varargin(i+1));
                else
                    plot(x,y);
                    hh=gce();
                    hh.children.foreground=CIindex; 
                end
                //______________.siso arrays.__________________________//
            elseif(or(typeof(varargin(i))==['rational','state-space']) & (typeof(varargin(i+1))=='boolean' )) then
                if(varargin(i+1)<>%T) then
                    error("Sigma: Siso array must be followed with boolean input %T");
                end
                xx=size(varargin(i),'r');
                yy=size(varargin(i),'c');
                zz=size(varargin(i),3);
                for ii=1:xx
                    for j=1:yy
                        for k=1:zz
                            y=svplot(varargin(i)(ii,j,k),x);
                            y=20*log10(y);
                            if (i<>rhs-1 & typeof(varargin(i+2))=='string') then
                                plot(x,y,varargin(i+2));
                            else
                                plot(x,y);
                                hh=gce();
                                hh.children.foreground=CIindex; 
                            end
                        end
                    end
                end
                
            end


        end
    end 
  xtitle("singular values","Frequency(rad/s)","Singular Values(dB)");
  H=gca();
    
    
endfunction
