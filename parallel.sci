function [varargout]=parallel(varargin)
    //function parallel
    //Parallel connection of two input/output models.
    //            |--------|
    //            |        |--------->z1
    // v1-------->|  SYS1  |
    //     u1 +-->|        |----+ y1
    //        |   |--------|    |
    //        |                 o------>y
    //  u---->+                 |
    //        |   |--------|    |
    //     u2 +-->|        |----+ y2
    //            |  SYS2  |
    //  v2------->|        |---------->z2
    //            |--------|
    //
    //SYS=PARALLEL(SYS1,SYS2) connects the input/output models SYS1 and SYS2 in 
    //    parallel.All the inputs SYS1 and SYS2 are connected correspondingly and corresponding 
    //    O/P are summed.If one of the system is SISO and other is MIMO then for every MIMO 
    //    corresponding system for i/p and o/p the SISO system is paralleled.
    //
    // SYS = PARALLEL(SYS1,SYS2,IN1,IN2,OUT1,OUT2) connects the input/output models  
    //      SYS1 and SYS2 in parallel. The inputs specified by IN1 and IN2 are connected 
    //      and the outputs specified by OUT1 and OUT2 are summed.The vectors IN1 and IN2 contain 
    //      indices into the input vectors of M1 and M2, respectively, and define 
    //      the input channels u1 and u2 in the diagram. Similarly, the vectors 
    //      OUT1 and OUT2 contain indexes into the outputs of M1 and M2. 
    //
    //[A,B,C,D] = PARALLEL(A1,B1,C1,D1,A2,B2,C2,D2)  produces a state-
    //             space system consisting of the parallel connection of systems 1 
    //             and 2 that connects all the inputs together and sums all the 
    //             outputs of the two systems,  Y = Y1 + Y2 (o/p of SYS1(A1,B1,C1,D1) and SYS2(A2,B2,C2,D2)).
    //
    //[A,B,C,D] = PARALLEL(A1,B1,C1,D1,A2,B2,C2,D2,IN1,IN2,OUT1,OUT2) 
 //                 connects the two systems in parallel such that the inputs 
//                 specified by IN1 and IN2 are connected and the outputs specified
//                 by OUT1 and OUT2 are summed. The vectors IN1 and IN2 contain 
//                 indexes into the input vectors of system 1 and system 2, 
//                 respectively.  Similarly, the vectors OUT1 and OUT2 contain 
//                 indexes into the outputs of the systems.  The parallel connection
//                 is performed by appending the two systems, summing the specified
//                 inputs and outputs, and removing the, now redundant, inputs and 
//                 outputs of system 2.
//
//  [NUM,DEN] = PARALLEL(NUM1,DEN1,NUM2,DEN2) produces a parallel 
//              connection of the two transfer function systems.
//
//  SYS=PARALLEL(SISOarray,'v',SYS2)
//  SYS=PARALLEL(SYS2,SISOarray,'v')
//                     connects the input/outputs of every element in siso array with SYS2 in 
//    parallel.   
//  
//  SYS=PARALLEL(SISOarray1,SISOarray2)
//              connects all the tf of first SISO array to second SISO array element wise in parallel.
//Examples:-
// s=%s;sys1=syslin('c',1/(500*s^2));sys2=syslin('c',(s+1)/(s+2));
//sys=parallel(sys1,sys2);              
//a1=[1 2;3 4];b1=[2 3;4 5];c1=[3 4;5 6];d1=[4 5;6 7];a2=[1 4;3 4];b2=[2 6;4 5];c2=[7 4;5 6];d2=[4 9;6 7];
//aa=syslin('c',a1,b1,c1,d1);bb=syslin('c',a2,b2,c2,d2);
// sys=parallel(aa,bb,1,2,2,1);
//[a b c d]=parallel(a1,b1,c1,d1,a2,b2,c2,d2,1,2,2,1);
//
//Author:-Paresh Yeole emailid:-yeoleparesh@students.vnit.ac.in
 
    [lhs,rhs]=argn(0);
    ni=length(varargin);
    nd=length(varargout);
    kk=1;
    
    select rhs
    //-----------------------------------------------------two systems case without vectors ---------------------------------------------//
     case 2 then
         if (and(typeof(varargin(1))<>['rational','state-space']) & and(typeof(varargin(2))<>['rational','state-space'])) then
           error("parallel:Wrong type of input arguments for given no. of input arguments");
         end
         for i=1:2
         if (varargin(i).dt=='d') then
            varargin(i).dt=1;
         end
         end
         if ((varargin(1).dt)<>(varargin(2).dt)) then
             error("parallel:sampling time must agree");
         end
         if(size(varargin(1))<>[1 1] & size(varargin(2))<>[1 1]) then
            if(size(varargin(1))<>size(varargin(2))) then
                 error("Incompatible sizes of the systems");
            else
                varargout(1)=varargin(1)+varargin(2);
            end
            
        elseif(size(varargin(1))==[1 1] & size(varargin(2))==[1 1]) then 
            varargout(1)=varargin(1)+varargin(2);
        elseif((size(varargin(1))==[1 1] | size(varargin(2))==[1 1]) & (typeof(varargin(1))=='state-space' | typeof(varargin(2))=='state-space')) then
            //disp("remain to code");
           [ny1,nu1] = size(varargin(1));
           [ny2,nu2] = size(varargin(2));
           if(size(varargin(1))==[1 1]) then
               if(typeof(varargin(1))=='state-space') then
                   sys=(varargin(1));
               else
                   sys=minss(tf2ss(varargin(1)));
               end
               
                a=sysdiag(sys.a,varargin(2).a);
                sys.b=sys.b*ones(size(sys.b,'c'),size(varargin(2).b,'c'));
                b=cat(1,sys.b,varargin(2).b);
                sys.c=ones(size(varargin(2).c,'c'),size(sys.c,'r'))*sys.c;
                c=cat(2,sys.c,varargin(2).c);
                d=varargin(2).d
            elseif(size(varargin(2))==[1 1]) then 
                if(typeof(varargin(2))=='state-space') then
                   sys=(varargin(2));
                else
                   sys=minss(tf2ss(varargin(2)));
               end
                a=sysdiag(varargin(1).a,sys.a);
                sys.b=sys.b*ones(size(sys.b,'c'),size(varargin(1).b,'c'));
                b=cat(1,varargin(1).b,sys.b);
//                disp(b);
                sys.c=ones(size(varargin(1).c,'c'),size(sys.c,'r'))*sys.c;
                c=cat(2,varargin(1).c,sys.c);
                d=varargin(1).d;
           end
         varargout(1)=syslin(varargin(1).dt,a,b,c,d);
        else 
            varargout(1)=varargin(1)+varargin(2);
         end
         

         //-----------------------------------------------two systems case with polynomials case---------------------------------//
    case 4 then
     num1 = coeff(varargin(1));
     den1 =coeff(varargin(2)); 
     num2 = coeff(varargin(3));
     den2 = coeff(varargin(4));
   if((length(num1)<=length(den1)) & (length(num2)<=length(den2)) ) then
     [nn,mn] = size(num1);
     for k=1:nn
       a(k,:) = conv(num1(k,:),den2) + conv(num2(k,:),den1);
       b = conv(den1,den2);
     end
     varargout(1)=a;
     varargout(2)=b;
   else
      error("parallel:input transfer function must be proper");
     end
    //----------------------------------------------two systems case with vector i/ps and o/ps--------------------------------------------//
    case 6 then
     /////case-parallel(sys1,sys2,inp1,inp2,out1,out2)
         if (and(typeof(varargin(1))<>['rational','state-space']) & and(typeof(varargin(2))<>['rational','state-space'])) then
           error("parallel:Wrong type of input arguments for given no. of input arguments");
         end
         for i=1:2
             if (varargin(i).dt=='d') then
                 varargin(i).dt=1;
             end
         end
         if ((varargin(1).dt)<>(varargin(2).dt)) then
             error("parallel:sampling time must agree");
         end
        [e f g h]=varargin(3:6);
         [ny1,nu1] = size(varargin(1));
        [ny2,nu2] = size(varargin(2));
        //State space systems with selection vectors
        if(e>ny1 | f>ny2) then
        error("parallel:specified inputs are out of range");
        
        end
        if(g>nu1 | h>nu2) then
        error("parallel:specified outputs are out of range");
        
        end
           inputs1 = e;      outputs1 = g;
           inputs2 = f+nu1;  outputs2 = h+ny1;
       
//        [a b c d]=abcd(varargin(1));
//        [a1 b1 c1 d1]=abcd(varargin(2));
//        
        //Check sizes
        if (length(inputs1)<>length(inputs2))  then
           error('Input sizes don''t match.');
        elseif (length(outputs1)<>length(outputs2)) then
           error('Output sizes don''t match.');
        end
         
        if(typeof(varargin(1))=='state-space' | typeof(varargin(2))=='state-space') then
                    if((typeof(varargin(1))<>'state-space') & (size(varargin(1),'*')<>1 | size(varargin(2),'*')<>1 ) ) then
                        varargin(1)=minss(tf2ss(varargin(1)));
                    elseif (typeof(varargin(2))<>'state-space' & (size(varargin(1),'*')<>1 | size(varargin(2),'*')<>1) ) then
                        varargin(2)=minss(tf2ss(varargin(2)));
                    end
                    sys=sysdiag(varargin(1),varargin(2));

            a=sys.a;b=sys.b;c=sys.c;d=sys.d;
            b(:,inputs1)=b(:,inputs1)+b(:,inputs2);
            d(:,inputs1)=d(:,inputs1)+d(:,inputs2);
            c(outputs1,:)=c(outputs1,:)+c(outputs2,:);
            d(outputs1,:)=d(outputs1,:)+d(outputs2,:);
            b(:,inputs2) = [];     d(:,inputs2) =  [];
            c(outputs2,:) = [];    d(outputs2,:) = [];
            varargout(1)=syslin('c',a,b,c,d);
        else
            sys=sysdiag(varargin(1),varargin(2));
            sys(:,inputs1)=sys(:,inputs1)+sys(:,inputs2);
            sys(outputs1,:)=sys(outputs1,:)+sys(outputs2,:);
            sys(:,inputs2) =  []; sys(outputs2,:) = [];
            varargout(1)=sys;
         end
         //-----------------------------------------------------------state-space case of 8------------------------------------------------------------------//
      case 8 then
         //disp("hi")
        ////state space systems  
         [a b c d]=varargin(1:4);
         [a1 b1 c1 d1]=varargin(5:8);
        if(size(b)<>[mtlb_size(a,1),mtlb_size(d,2)] | size(c)<>[mtlb_size(d,1),mtlb_size(a,2)] | size(d)<>[mtlb_size(c,1),mtlb_size(b,2)] | issquare(a)==%F) then
               error('wrong size of matrices given for state space model');
        end
        if(size(b1)<>[mtlb_size(a1,1),mtlb_size(d1,2)] | size(c1)<>[mtlb_size(d1,1),mtlb_size(a1,2)] | size(d1)<>[mtlb_size(c1,1),mtlb_size(b1,2)] | issquare(a1)==%F) then
               error('wrong size of matrices given for state space model');
        end
        [ny1,nu1] = size(d);
        [ny2,nu2] = size(d1);
         // State space systems w/o selection vectors
         inputs1 = [1:nu1];     outputs1 = [1:ny1];
         inputs2 = [1:nu2]; outputs2 = [1:ny2]; 
       
        //Check sizes
        if (length(inputs1)<>length(inputs2))  then
           error('Input sizes don''t match.');
        elseif (length(outputs1)<>length(outputs2)) then
           error('Output sizes don''t match.');
        end
        sys1=syslin('c',a,b,c,d);
        sys2=syslin('c',a1,b1,c1,d1);
        sys=sys1+sys2;
        varargout(1)=sys.a;
        varargout(2)=sys.b;
        varargout(3)=sys.c;
        varargout(4)=sys.d;
        //-------------------------------------------------------Case of 12 - state-space matrices--------------------------------------------------------------------------------//
    case 12 then
        [a b c d]=varargin(1:4);
         [a1 b1 c1 d1]=varargin(5:8);
        if(size(b)<>[mtlb_size(a,1),mtlb_size(d,2)] | size(c)<>[mtlb_size(d,1),mtlb_size(a,2)] | size(d)<>[mtlb_size(c,1),mtlb_size(b,2)] | issquare(a)==%F) then
               error('wrong size of matrices given for state space model');
        end
        if(size(b1)<>[mtlb_size(a1,1),mtlb_size(d1,2)] | size(c1)<>[mtlb_size(d1,1),mtlb_size(a1,2)] | size(d1)<>[mtlb_size(c1,1),mtlb_size(b1,2)] | issquare(a1)==%F) then
               error('wrong size of matrices given for state space model');
        end
        [ny1,nu1] = size(d);
        [ny2,nu2] = size(d1);
        [e f g h]=varargin(9:12);
        if(e>ny1 | f>ny2) then
        error("parallel:specified inputs are out of range");
        
        end
        if(g>nu1 | h>nu2) then
        error("parallel:specified outputs are out of range");
        
        end
        //State space systems with selection vectors
           inputs1 = e;      outputs1 = g;
           inputs2 = f+nu1;  outputs2 = h+ny1;
        
        //Check sizes
        if (length(inputs1)<>length(inputs2))  then
           error('Input sizes don''t match.');
        elseif (length(outputs1)<>length(outputs2)) then
           error('Output sizes don''t match.');
        end
        sys1=syslin('c',a,b,c,d);
        sys2=syslin('c',a1,b1,c1,d1);
        //sys=sys1+sys2;
        
        sys=sysdiag(sys1,sys2);a=sys.a;b=sys.b;c=sys.c;d=sys.d;
        b(:,inputs1)=b(:,inputs1)+b(:,inputs2);
        d(:,inputs1)=d(:,inputs1)+d(:,inputs2);
        c(outputs1,:)=c(outputs1,:)+c(outputs2,:);
        d(outputs1,:)=d(outputs1,:)+d(outputs2,:);
         b(:,inputs2) = [];     d(:,inputs2) =  [];
         c(outputs2,:) = [];    d(outputs2,:) = [];
        varargout(1)=a;
        varargout(2)=b;
        varargout(3)=c;
        varargout(4)=d;
    case 3 then
         if (and(typeof(varargin(1))<>['rational','state-space']) & and(typeof(varargin(2))<>['rational','state-space','string']) & and(typeof(varargin(3))<>['rational','state-space','string'])) then
           error("parallel:Wrong type of input arguments and string v  is expected to pass  with SISO array with MIMO/SISO ");
         end
         
//         if(typeof(varargin(2))<>'rational') then
//         error("parallel:to pass SISO array with MIMO SISO array must be the 2nd i/p argument");
//         end
         
         for i=1:3
             if(typeof(varargin(i))=='string') then
                 sisoindex=i;
                 continue;
             end
             
            if (varargin(i).dt=='d') then
                varargin(i).dt=1;
            end
         end
         if(varargin(sisoindex)<>'v') then
             error("parallel:string v is expected as after SISO array with MIMO/SISO systems");
         end
         sisoarray=varargin(sisoindex-1);
         if(sisoindex==2) then
             
         if ((varargin(1).dt)<>(varargin(3).dt)) then
             error("parallel:sampling time must agree");
         end
         mimo=varargin(3);
         if(typeof(varargin(3))=='rational') then
            mimo=minss(tf2ss(varargin(3)));
         end
         elseif(sisoindex==3) then
             
         if ((varargin(1).dt)<>(varargin(2).dt)) then
             error("parallel:sampling time must agree");
         end
         mimo=varargin(1);
         if(typeof(varargin(1))=='rational') then
            mimo=minss(tf2ss(varargin(1)));
         end
         end
         
         
        for i=1:size(sisoarray,'r')
            for j=1:size(sisoarray,'c')
                for k=1:size(sisoarray,3)
                sys=minss(tf2ss(sisoarray(i,j,k)));
                a=sysdiag(mimo.a,sys.a);
                sys.b=sys.b*ones(size(sys.b,'c'),size(mimo,'c'));
                b=cat(1,mimo.b,sys.b);
                sys.c=ones(size(mimo,'c'),size(sys.c,'r'))*sys.c;
                c=cat(2,mimo.c,sys.c);
                d=sys.d+mimo.d;
                t=syslin(varargin(1).dt,a,b,c,d);
                disp(msprintf(gettext("OUTPUT--%d*%d*%d"),i,j,k))
                
                disp(t);
                //kk=kk+1;
                end
            end
        end
        varargout(1)=(msprintf(gettext("%d*%d*%d state-space array"),size(sisoarray,'r'),size(sisoarray,'c'),size(sisoarray,3)));
   else
       error("Wrong No. of i/p arguments");
   end;

endfunction
//
//_________________________________________________________
