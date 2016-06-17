function [o]=impulseplot(varargin)

//Author:- Paresh Yeole emailid:-yeoleparesh@students.vnit.ac.in
//Gives the impulse response of continuos and discrete SISO as well as MIMO  systems  
//impulseplot  Plot impulse response of linear systems.
//impulseplot, an extension of IMPULSE, provides a command line interface for customizing the plot appearance.
//
//calling sequence:-
//impulse(sys)
//impulse(poly1,poly2)
//impulse(sys,Tfinal)
//impulse(sys,Tvector)
//impulse(sys1,sys2,...,T)
//impulse(sys1,'r',sys2,'y--',sys3,'gx',..)
//
//parameters:-
//sys:- sys can be SISO array,MIMO system or SISO either discrete or continuous
//poly1:-numerator of the system
//poly2:-denominator of the system
//Tfinal:-time upto which the response is to be calculated/plotted.
//Tvector:-time vector through which the response is to be plotted.For Discrete time, the sampling time must match the spacing of the given time vector.
//o -handle of the figure plotted to customize the plot appearances
//
//
// Note: In discrete time, IMPULSE computes the response to a unit-area 
//  pulse of length Ts and height 1/Ts where Ts is the sample time. This
//  pulse approaches the continuous-time Dirac impulse delta(t) as Ts goes
//  to zero.
////Examples:-
//s=poly(0,'s');
//sys=syslin('c',(s+3)/(s^3+4*s+2));
//impulseplot(sys)   plots the impulse response of sys
//sys1=ssrand(2,3,4);
//impulseplot(sys1)
//impulseplot(sys,sys1)
//impulseplot(sys,'--r',sys1,'gx')
//aa=pid(rand(2,3,4),2,3,4);
//impulseplot(aa,%T)
  
    [lhs,rhs]=argn(0);
  
    ni=length(varargin);
    flag=0;
       
    
    
    if rhs == 0 | (rhs == 1 & typeof(varargin($)) <> ['state-space', 'rational'] | (rhs == 1 & size(varargin($)) == [0 0])) then
        error(msprintf(gettext("%s: Wrong type for input argument \n#%d: State-space or transfer function of linear system expected or two polynomials expected.\n"),"impulse",1))
    end
    if typeof(varargin(1)) <> ['state-space','rational','polynomial'] then
        error(msprintf(gettext("%s: Wrong type for first input argument\n#%d: State-space or transfer function expected.\n"),"impulse",1))
    end
    
    
    t=0:0.1:100;
    
    
   if rhs > 1 then
        if typeof(varargin($)) == 'constant' then
            if size(varargin($)) == [1 1] then
                if varargin($) <= 0 then
                    error(msprintf(gettext("%s: The final time value must be a positive real number.\n"),"impulse"))
                end
                tFinal = varargin($)
                t =0:0.01:tFinal; 
                flag=1;
            elseif isequal(size(varargin($)),[1 1]) == %f then
                // finding that the time vector has positive time value
                if(and(diff(varargin($)<0))) then  //if time vector is not monotonically increasing
                    error(msprintf("impulse: the time vector must be real, finite, and must contain monotonically increasing and evenly spaced time samples."));
                end
                flag=1;
                tempTimeIndex = find(varargin($) >= 0) //finding the positive no.'s index in the vector
                if isequal(size(varargin($)),size(tempTimeIndex)) == %t then
                    temp=varargin($);
                    if ((typeof(varargin(1))<>['polynomial']) & (varargin(1).dt<>'c')) then
                        if(varargin(1).dt=='d') then
                            dt=1;
                        else
                            dt=varargin(1).dt;
                            end
                            if(dt<>(temp(2)-temp(1))) then
                            error(msprintf("impulse: for discrete time system the vector difference must be equal to the sample time of the system"));
                            end
                            
                        end
                                          t = varargin($)
                else
                    tempTime=varargin($);
                    tempTime = tempTime(tempTimeIndex(1):tempTimeIndex($))
                    t = tempTime
                end
            end
        end
    end
    
    
    ///////////////////////////////////////////////////
    //getting the sublot credentials
    for i=1:rhs
        ttf=varargin(i);
        if typeof(ttf)=='state-space' then
                ttf=ss2tf(ttf);
            end
            if typeof(ttf)=='rational' then
        xx1(i)=size(ttf,'r');
        yy1(i)=size(ttf,'c');
        end
    end
    
    /////////////////////////////////printf("\n%d,%d",max(xx1),max(yy1));
    
    
    
     CIindex=1;
   
        for i=1:rhs
            CIindex=CIindex+1;
            I=1;
            if(or(typeof(varargin(i))==['rational','state-space']) & ((varargin(i).dt)<>'c')) then
            if((typeof(varargin($))=='constant') & (size(varargin($)) <> [1 1]) & ((varargin(i).dt)<>(t(2)-t(1)))) then
            error(msprintf("sampling time of the given system must match the step of the vector"));
        elseif((typeof(varargin($))=='constant') & (size(varargin($)) == [1 1])) then
            
            temptime=t;
            if(varargin(i).dt=='d') then
                dt=1
            else
                dt=varargin(i).dt
            end
            
            t=0:varargin(i).dt:varargin($);
        else
            //temptime=t;
            if(varargin(i).dt=='d') then
                dt=1;
            else
               
                dt=varargin(i).dt;
            end
            
            t=0:varargin(i).dt:t($);
       
            
            end
            
            end
            
            if typeof(varargin(i))=='state-space' then
                varargin(i)=ss2tf(varargin(i));
            end
            
            /////////////////////////SISO system//////////////////////////////
        if((typeof(varargin(i))=='rational') & (size((varargin(i)),'*')==1)) then
        //sysI=sysI+1;
        
         if flag<>1 then
             pp=pole(varargin(i))
             pp=cell2mat(pp);
              ppr=real(pp)
               if(varargin(i).dt=='c') then
                y=csim('impuls',t,varargin(i));
            else
                y=flts(eye(1,length(t)),varargin(i));
            end
             if or(ppr > 0) then
                 //ch=find(y>=10^12 | y<=-10^12);
                 dompol=max(ppr);
                 tfinal=25/(dompol*log10(%e));
                  if((varargin(i).dt)<>'c') then
                     //t=0:varargin(i).dt:t(ch(1));
                     t=0:varargin(i).dt:tfinal;
                  else
                     //t=0:0.1:t(ch(1))
                     t=0:0.1:tfinal;

                  end
             elseif and(ppr<=0) then
                
                //y=csim('impuls',t,varargin(i));
                 for iii=length(t):-1:1
                  if(y(iii)<-0.002 | y(iii)>0.002) then
                      break; 
                  end;
                  
              end
              if((varargin(i).dt)<>'c') then
                  t=0:varargin(i).dt:(iii-1)*0.1;
                  else
                      
                   t=0:0.1:(iii-1)*0.1;
                  end
                     
                end
             end
             
         if i<>rhs & typeof(varargin(i+1))=='string'  then
           if (varargin(i).dt=='c') then
           G=csim('impuls',t,varargin(i));
       else
          G=flts(eye(1,length(t)),varargin(i));
       end
       subplot(max(xx1),max(yy1),I);
       plot(t,G,varargin(i+1));
       if (varargin(i).dt<>'c') then
          hh=gce();
         hh.children.polyline_style=2;  
         end
        else
            if (varargin(i).dt=='c') then
           G=csim('impuls',t,varargin(i));
       else
          G=flts(eye(1,length(t)),varargin(i));
       end
       subplot(max(xx1),max(yy1),I);
       plot(t,G);
       hh=gce();
       hh.children.foreground=CIindex;
       if (varargin(i).dt<>'c') then
       
         hh.children.polyline_style=2;  
       end
      
        end
       
    
        ////////////////SISO array///////////////////////////

     elseif typeof(varargin(i))=='rational' & size(varargin(i),'*')<>1 &  i<>rhs & rhs<>1 & typeof(varargin(i+1))=='boolean' then
         if(varargin(i+1)<>%T   ) then
             error(msprintf("impulse:wrong input arguments"));
         end
         
        xx=size(varargin(i),'r');
        yy=size(varargin(i),'c');
        zz=size(varargin(i),3);
        tt=varargin(i);
        I=1;
        if flag<>1 then
            for ii=1:xx
                 for jj=1:yy
                     for kk=1:zz 
             pp=pole(tt(ii,jj))
             pp=cell2mat(pp);
              ppr=real(pp)
                 if(varargin(i).dt=='c') then
                y=csim('impuls',t,tt(ii,jj,kk));
            else
                y=flts(eye(1,length(t)),tt(ii,jj,kk));
            end
              if or(ppr > 0) then
                 //ch=find(y>10^30);
                  dompol=max(ppr);
                 tfinal=25/(dompol*log10(%e));
                 temp(ii,jj)=tfinal;
                  //temp=t(ch(1));
             elseif and(ppr<=0) then
                for iii=length(t):-1:1
                  if(y(iii)<-0.002 | y(iii)>0.002) then
                      break; 
                  end;
                  
                  end
                    temp(ii,jj)=(iii-1)*0.1;
                    end 
                end
                end
            end
            if((varargin(i).dt)<>'c') then
                  t=0:varargin(i).dt:max(temp);
                  else
                      
                   t=0:0.1:max(temp);
                  end
            //t=0:0.1:max(temp);
             end
        if i<rhs-1 & typeof(varargin(i+2))=='string'  then
        for ii=1:xx
          for jj=1:yy
              for kk=1:zz
              if(varargin(i).dt=='c') then
          G=csim('impuls',t,tt(ii,jj,kk));
      else
          G=flts(eye(1,length(t)),tt(ii,jj,kk));
          end
         title(msprintf(gettext("input(1)-output(1)")));
         plot(t,G,varargin(i+2));
         hh=gce();
         if(varargin(i).dt<>'c') then
         
         hh.children.polyline_style=2;
         end
     
         //I=I+1;
          end
          end
        end
        
    else
        for ii=1:xx
          for jj=1:yy
              for kk=1:zz
              if(varargin(i).dt=='c') then
          G=csim('impuls',t,tt(ii,jj,kk));
      else
          G=flts(eye(1,length(t)),tt(ii,jj,kk))
          end
         title(msprintf(gettext("input(1)-output(1)")));
        // ylabel("output(%d)",jj);
         plot(t,G);
         hh=gce();
         hh.children.foreground=CIindex;
         if(varargin(i).dt<>'c') then
         
         hh.children.polyline_style=2;
        end
         //I=I+1;
          end
          end
        end
        end
        
        
   /////////////////////////////MIMO-system////////////////////////
    elseif ((typeof(varargin(i))=='rational') & (size(varargin(i),'*')<>1)) then
        xx=size(varargin(i),'r');
        yy=size(varargin(i),'c');
        tt=varargin(i);
        I=1;
        if flag<>1 then
            for ii=1:xx
                 for jj=1:yy
             pp=pole(tt(ii,jj))
             pp=cell2mat(pp);
              ppr=real(pp)
               if(varargin(i).dt=='c') then
                y=csim('impuls',t,tt(ii,jj));
            else
                y=flts(eye(1,length(t)),tt(ii,jj));
            end
              if or(ppr > 0) then
                  //ch=find(y>=10^30);
                  //temp=t(ch(1)-1);
                  dompol=max(ppr);
                 tfinal=25/(dompol*log10(%e));
                 temp(ii,jj)=tfinal;
             elseif and(ppr<=0) then
                for iii=length(t):-1:1
                  if(y(iii)<-0.002 | y(iii)>0.002) then
                      break; 
                  end;
//                  
                  end
                    temp(ii,jj)=(iii-1)*0.1;
                    end 
                end
            end
            if((varargin(i).dt)<>'c') then
                  t=0:varargin(i).dt:max(temp);
                  else
                      
                   t=0:0.1:max(temp);
                  end
             end
        if i<>rhs & typeof(varargin(i+1))=='string'  then
        for ii=1:xx
          for jj=1:yy
              if(varargin(i).dt=='c') then
          G=csim('impuls',t,tt(ii,jj));
      else
          G=flts(eye(1,length(t)),tt(ii,jj));
          end
         subplot(max(xx1),max(yy1),I)
         title(msprintf(gettext("input(%d)-output(%d)"),ii,jj));
         plot(t,G,varargin(i+1));
        if(varargin(i).dt<>'c') then
        hh=gce();
         hh.children.polyline_style=2;
         end
         I=I+1;
          end
          
        end
    else
        for ii=1:xx
          for jj=1:yy
              if(varargin(i).dt=='c') then
          G=csim('impuls',t,tt(ii,jj));
      else
          G=flts(eye(1,length(t)),tt(ii,jj))
          end
         subplot(max(xx1),max(yy1),I)
         title(msprintf(gettext("input(%d)-output(%d)"),ii,jj));
         plot(t,G);
         hh=gce();
         hh.children.foreground=CIindex;
         if(varargin(i).dt<>'c') then
         
         hh.children.polyline_style=2;
         end
         I=I+1;
          end
           end
        end
    end
    end
       
   
    ///////////////////////////////////////////////////////////////////////////////////////
    if(typeof(varargin(1))=='polynomial') then
   if(typeof(varargin(2))=='polynomial') then
        if flag<>1 then
             pp=pole(varargin(1)/varargin(2))
             pp=cell2mat(pp);
              ppr=real(pp)
             if or(ppr >= 0) then
                  dompol=max(ppr);
                 tfinal=25/(dompol*log10(%e));
                 t=0:0.1:tfinal;
             elseif and(ppr<0) then
                y=csim('impuls',t,varargin(1)/varargin(2));
                 for iii=length(t):-1:1
                  if(y(iii)<-0.002 | y(iii)>0.002) then
                      break; 
                  end;
//                  if(y(iii)<=10^4) then
//                      break;
//                  end;
                  
                  end
                    t=0:0.1:(iii-1)*0.1; 
                end
             end
   G=csim('impuls',t,varargin(1)/varargin(2));
   
        if(typeof(varargin(3))=='string') then
             plot(t,G,varargin(3));
             else
             
         plot(t,G);
  end 
else
   error('not enough arguments'); 
end
end
    /////////////////////////////////////////////////////
    
    o=gcf();
    o.figure_name='IMPULSE-PLOT'
    
    
    
endfunction
