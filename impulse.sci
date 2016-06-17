function[varargout]=impulse(varargin)
    //Author:- Paresh Yeole emailid:-yeoleparesh@students.vnit.ac.in
//Gives the impulse response of continuos and discrete SISO as well as MIMO  systems  
//
//calling sequence:-
//impulse(sys)
//impulse(poly1,poly2)
//impulse(sys,Tfinal)
//impulse(sys,Tvector)
//impulse(SISOarray,%T)
//impulse(sys1,sys2,...,T)
//impulse(sys1,'r',sys2,'y--',sys3,'gx',..)
//[Y,t]=impulse(sys)
//For state-space models, 
//[Y,T,X] = impulse(SYS, ...) 
//Response uncertainty computation: 
//[Y,T,X,YSD] = impulse(SYS)
//
//parameters:-
//sys:- sys can be SISO array,MIMO system or SISO either discrete or continuous
//poly1:-numerator of the system
//poly2:-denominator of the system
//Tfinal:-time upto which the response is to be calculated/plotted.
//Tvector:-time vector through which the response is to be plotted.For Discrete time, the sampling time must match the spacing of the given time vector.
//Y:-vector of the impulse response.
//T:-time vector through which response is calculate.
//X:-state-Trajectory response vector.
//YSD:- standard deviation YSD of the   response Y of an identified system SYS. YSD is empty if SYS does not  contain parameter covariance information.
//
// Note: In discrete time, IMPULSE computes the response to a unit-area 
//  pulse of length Ts and height 1/Ts where Ts is the sample time. This
//  pulse approaches the continuous-time Dirac impulse delta(t) as Ts goes
//  to zero.
//Examples:-
//s=poly(0,'s');
//sys=syslin('c',(s+3)/(s^3+4*s+2));
//impulse(sys)   plots the impulse response of sys
//[y t]=impulse(sys)   gives y as a impulse response matrix and t as time vector through which the response is calculated
//sys1=ssrand(2,3,4);
//impulse(sys1)
//[y t x ysd]=impulse(sys1);
//impulse(sys,sys1)
//impulse(sys,'--r',sys1,'gx')
//aa=pid(rand(2,3,4),2,3,4);
//impulse(aa,%T)


    [lhs,rhs]=argn(0);
    
    nd=length(varargout);
    ni=length(varargin);
    flag=0;
    if(nd>1)
        if(ni>=2) then
            if(ni==2) then
            if or(typeof(varargin(2)) == ['rational', 'state-space','string']) then
            error(msprintf("The impulse command operates on a single model when used with output arguments."));
            end
        else
            error(msprintf("The impulse command operates on a single model when used with output arguments."));
            end
            
        end
        
    
    end
    
    
    
    
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
    if lhs==1 then
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
                 //disp(ch);
                 dompol=max(ppr);
                 tfinal=25/(dompol*log10(%e));
                  if((varargin(i).dt)<>'c') then
                    // t=0:varargin(i).dt:t(ch(1));
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
             varargout(1)="IMPULSE_PLOT";
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
     elseif typeof(varargin(i))=='rational' & size(varargin(i),'*')<>1 & i<>rhs & rhs<>1 & typeof(varargin(i+1))=='boolean' then
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
                  dompol=max(ppr);
                 tfinal=25/(dompol*log10(%e));
                 temp(ii,jj)=tfinal;
//                 ch=find(y>10^8 | y<-10^8);
//                  //temp=100;
//                  temp=t(ch(1));
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
        
        varargout(1)="IMPULSE_PLOT";
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
                  //ch=find(y>=10^12 | y<=-10^12);
                  //temp=t(ch(1));
                  //disp(temp);
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
    varargout(1)="IMPULSE_PLOT";
        end
    h=gcf();
    h.figure_name= "IMPULSE-PLOT"; 
    
    ///////////////////more than 1 o/p arguments/////////////////////////////////////////////////////////
elseif lhs>1 then
        
   //varargout(1)=(csim('impuls',t,varargin(1)))';
    xx=size(varargin(1),'r');
        yy=size(varargin(1),'c');
        op1=cell(xx,yy);
        op2=cell(xx,yy);
        tt=varargin(1);
         for ii=1:xx
                 for jj=1:yy
                     if(varargin(1).dt=='c') then
      [G x1]=csim('impuls',t,tt(ii,jj));
  else
      if(typeof(varargin(1))=='state-space') then
      [G x1]=flts(eye(1,length(t)),tt(ii,jj));
  else
      [G]=flts(eye(1,length(t)),tt(ii,jj));
  end
  
      end
      op1(ii,jj).entries=G';
      if(typeof(varargin(1))=='state-space') then
      op2(ii,jj).entries=x1';
      end

 
end
end
varargout(1)=cell2mat(op1);
varargout(2)=t';

if(typeof(varargin(1))=='state-space') then
           varargout(3)=cell2mat(op2);
            varargout(4)=stdev(varargout(1));
       else
            varargout(3)=[];
            varargout(4)=[];
        end
 
    end
    
    
//////////////////////////if both initial arguments are polynomials////////////////    
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
                  
                  end
                    t=0:0.1:(iii-1)*0.1; 
                end
             end
   G=csim('impuls',t,varargin(1)/varargin(2));
   varargout(1)=(G)';
         if(typeof(varargin(3))=='string') then
             plot(t,G,varargin(3));
             else
             
         plot(t,G);
  end 
else
   error('not enough arguments'); 
end
end

    
    
    

endfunction
