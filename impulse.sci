//Author:- Paresh Yeole emailid:-yeoleparesh@students.vnit.ac.in
//











function[varargout]=impulse(varargin)
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
    
    printf("\n%d,%d",max(xx1),max(yy1));
    
    
    if lhs==1 then
        for i=1:rhs
            I=1;
            if(varargin(i).dt<>'c') then
            if((typeof(varargin($))=='constant') & (size(varargin($)) <> [1 1]) & varargin(i).dt<>(t(2)-t(1))) then
            error(msprintf("sampling time of the given system must match the step of the vector"));
        else
            temptime=t;
            if(varargin(i).dt=='d') then
                dt=1
            else
                dt=varargin(i).dt
            end
            
            t=0:varargin(i).dt:temptime($);
            end
            
            end
            
            if typeof(varargin(i))=='state-space' then
                varargin(i)=ss2tf(varargin(i));
            end
        if typeof(varargin(i))=='rational' & size(varargin(i),'*')==1 then
        //sysI=sysI+1;
        
         if flag<>1 then
             pp=pole(varargin(i))
             pp=cell2mat(pp);
              ppr=real(pp)
             if or(ppr >= 0) then
                  t=0:0.1:100
             elseif and(ppr<0) then
                 if(varargin(i).dt=='c') then
                y=csim('impuls',t,varargin(i));
            else
                y=flts(eye(1,length(t)),varargin(i));
            end
                //y=csim('impuls',t,varargin(i));
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
             if (i==1) then
                 if (varargin(i).dt=='c') then
                     varargout(1)=(csim('impuls',t,varargin(1)))';
                 else
                     varargout(1)=(flts(eye(1,length(t)),varargin(i)))';
                 end
            // varargout(1)=(csim('impuls',t,varargin(1)))';
         end
         if i<>rhs & typeof(varargin(i+1))=='string'  then
           if (varargin(i).dt=='c') then
           G=csim('impuls',t,varargin(i));
       else
          G=flts(eye(1,length(t)),varargin(i));
       end
       subplot(max(xx1),max(yy1),I);
         plot(t,G,varargin(i+1));
        else
            if (varargin(i).dt=='c') then
           G=csim('impuls',t,varargin(i));
       else
          G=flts(eye(1,length(t)),varargin(i));
       end
       subplot(max(xx1),max(yy1),I);
         plot(t,G);
        end
        
    
        
     elseif typeof(varargin(i))=='rational' & size(varargin(i),'*')<>1 & typeof(varargin(i+1))=='boolean' then
         if(varargin(i+1)<>%T   ) then
             error(msprintf("impulse:wrong input arguments"));
         end
         
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
              if or(ppr >= 0) then
                  temp=100;
             elseif and(ppr<0) then
             
                     if(varargin(i).dt=='c') then
                y=csim('impuls',t,tt(ii,jj));
            else
                y=flts(eye(1,length(t)),tt(ii,jj));
            end
            
                 for iii=length(t):-1:1
                  if(y(iii)<-0.002 | y(iii)>0.002) then
                      break; 
                  end;
//                  if(y(iii)<=10^4) then
//                      break;
//                  end;
//                  
                  end
                    temp(ii,jj)=(iii-1)*0.1;
                    end 
                end
            end
            t=0:0.1:max(temp);
             end
        if i<>rhs & typeof(varargin(i+2))=='string'  then
        for ii=1:xx
          for jj=1:yy
              if(varargin(i).dt=='c') then
          G=csim('impuls',t,tt(ii,jj));
      else
          G=flts(eye(1,length(t)),tt(ii,jj));
          end
         //subplot(max(xx1),max(yy1),I)
         title(msprintf(gettext("input(1)-output(1)")));
         plot(t,G,varargin(i+2));
         //I=I+1;
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
         //subplot(max(xx1),max(yy1),I)
         //vv=msprintf(gettext("input(%d)-output(%d)"),ii,jj);
         title(msprintf(gettext("input(1)-output(1)")));
        // ylabel("output(%d)",jj);
         
         plot(t,G);
         //I=I+1;
          end
          
        end
        end
        
        
   
    elseif typeof(varargin(i))=='rational' & size(varargin(i),'*')<>1 then
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
              if or(ppr >= 0) then
                  temp=100;
             elseif and(ppr<0) then
             
                     if(varargin(i).dt=='c') then
                y=csim('impuls',t,tt(ii,jj));
            else
                y=flts(eye(1,length(t)),tt(ii,jj));
            end
            
                 for iii=length(t):-1:1
                  if(y(iii)<-0.002 | y(iii)>0.002) then
                      break; 
                  end;
//                  if(y(iii)<=10^4) then
//                      break;
//                  end;
//                  
                  end
                    temp(ii,jj)=(iii-1)*0.1;
                    end 
                end
            end
            t=0:0.1:max(temp);
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
         //vv=msprintf(gettext("input(%d)-output(%d)"),ii,jj);
         title(msprintf(gettext("input(%d)-output(%d)"),ii,jj));
        // ylabel("output(%d)",jj);
         
         plot(t,G);
         I=I+1;
          end
           end
        end
        end
        end
   // end
    //lhs > 1
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
//      msprintf("%d %d %d",G);
 
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
    
    
    
//    pp=pole(sys)
//    ppr=real(pp)
//    if or(ppr>0) then
//        t=0:0.1:100
//    elseif and(ppr<=0) then
//        t=0:0.1:10^8;
//        l=1+0.2;
//        h=1-0.2;
//        u=ones(1,length(t));
//        y=csim('impuls',t,sys);
//        for iii=length(t):-1:1
//            if(y(iii)<l | y(iii)>h) break; end;
//            end
//       t=(iii-1)*0.1; 
//    end
    //G=csim('impuls',0:0.1:1000,sys);
    
    
    // if num,den is given then
///////////////////////////////////////////////////////////////////////////////////////////////////////
if(typeof(varargin(1))=='polynomial') then
   if(typeof(varargin(2))=='polynomial') then
        if flag<>1 then
             pp=pole(varargin(1)/varargin(2))
             pp=cell2mat(pp);
              ppr=real(pp)
             if or(ppr >= 0) then
                  t=0:0.1:100
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
   ////////////////////////////////////////////////////////////////////////////////////////
    
    
    
       ////////////////////////
       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
//    if rhs == 0 | (rhs == 1 & typeof(varargin($)) <> ['state-space', 'rational'] | (rhs == 1 & size(varargin($)) == [0 0])) then
//        error(msprintf(gettext("%s: Wrong type for input argument \n#%d: State-space or transfer function of linear system expected.\n"),"impulse",1))
//    end
//    if typeof(varargin(1)) <> ['state-space','rational'] then
//        error(msprintf(gettext("%s: Wrong type for first input argument\n#%d: State-space or transfer function expected.\n"),"impulse",1))
//    end
//    
//    //determining the time 
//    t = 0:0.1:1000-0.1; 
//    if rhs > 1 then
//        if typeof(varargin($)) == 'constant' then
//            if size(varargin($)) == [1 1] then
//                if varargin($) <= 0 then
//                    error(msprintf(gettext("%s: The final time value must be a positive real number.\n"),"impulse"))
//                end
//                tFinal = varargin($)
//                t =0:0.01:tFinal; 
//            elseif isequal(size(varargin($)),[1 1]) == %f then
//                // finding that the time vector has positive time value
//                if(and(diff(varargin($)<0))) then  //if time vector is not monotonically increasing
//                    error(msprintf("impulse: the time vector must be real, finite, and must contain monotonically increasing and evenly spaced time samples."));
//                end
//                
//                tempTimeIndex = find(varargin($) >= 0) //finding the positive no.'s index in the vector
//                if isequal(size(varargin($)),size(tempTimeIndex)) == %t then
//                    t = varargin($)
//                else
//                    tempTime=varargin($);
//                    tempTime = tempTime(tempTimeIndex(1):tempTimeIndex($))
//                    t = tempTime
//                end
//            end
//        end
//    end
//    
//    //for n systems and formatting
//    if rhs >= 1 then
//        if typeof(varargin($)) == 'constant' then
//            vararginIndex = rhs-1
//        else
//            vararginIndex = rhs
//        end
//    end
//    
//    DI=1;//DATA-index-no.
//    currentI = 1
//    previousI = 0
//    colorIndexNumber = 1
//    colorIndex = []
//    
//    for i=1:vararginIndex
//        printf("cheenu");
//        if typeof(varargin(i)) == 'state-space' | typeof(varargin(i)) == 'rational' then
//        D(DI,1)=i;
//        tempsize=size(varargin(i));
//        if typeof(varargin(ii)) == 'state-space' then
//                outputSize(DI,1) = tempsize(1,1) 
//                inputSize(DI,1) = tempsize(1,2)
//             elseif typeof(varargin(ii)) == 'rational' then
//                 if isequal(tempsize,[1 1]) == %f & length(tempsize) == 2 then
//                     outputSize(DI,1) = tempsize(1,1) 
//                     inputSize(DI,1) = tempsize(1,2)
//                 else
//                     outputSize(DI,1) = 1 
//                     inputSize(DI,1) = 1
//                 end
//        end
//        
//        if i==rhs then
//            tempi=rhs
//        else
//             tempi=i+1;
//        end
//        
//        if typeof(varargin(tempi)) == 'string' then
//        currentI = i+1
//                        if currentI == previousI+1 then
//                            error(msprintf(gettext("%s: Incorrect syntax.\n"),"stepplot"))
//                        end
//                    previousI = currentI
//                    colorIndex(DI,1) = ii+1
//                else
//                    colorIndex(DI,1) = 0
//                 end
//             
//         
//        end
//        DI = DI + 1
//    end
//    end
//    
//    if isequal(zeros(length(D),1),colorIndex) == %f then
//        [colorIndex newIndex] = gsort(colorIndex,'r')
//        dataIndex = dataIndex(newIndex(:,1),1)
//    end
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    //for variable outputs
//   
//    //t=0:0.1:10;
//    //p=diracdelta(t);
//    select argn(1)
//     case 1 then
//         varargout(1)=(csim('impuls',t,varargin(1)))';
//         G=csim('impuls',t,varargin(1));
//         plot(t,G);
// 
//    case 2 then
//        varargout(1)=(csim('impuls',t,varargin(1)))';
//        varargout(2)=t';
//    case 3 then
//        varargout(1)=(csim('impuls',t,varargin(1)))';
//       varargout(2)=t';
//        if(typeof(sys)=='state-space') then
//            varargout(3)=1;
//        else
//            varargout(3)=[];
//        end
//    case 4 then
//        varargout(1)=(csim('impuls',t,varargin(1)))';
//       varargout(2)=t';
//        if(typeof(sys)=='state-space') then
//            varargout(3)=1;
//            varargout(4)=stdev(varargout(1));
//        else
//            varargout(3)=[];
//            varargout(4)=[];
//        end
//        
//    end
    
   // end
   //G=csim('impuls',t,sys);
//diw=(((horner(sys.den,%i*2*%pi*(1/t)))))*(conj(horner(sys.den,%i*2*%pi*(1/t))));
//mag=sqrt(nn.*mm);
//plot(t,G);
endfunction
