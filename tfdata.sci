function[num,den,Ts]=tfdata(p)
    //author:- Paresh Yeole
//[num,den,Ts] = tfdata(sys) returns coefficients of the numerator(s) and denominator(s) of the transfer function for the TF, SS and the sampling time of the system.
//function -- tfdata(TRANSFER FUNCTION)
//sampling time of the system can also be find out.
//SISO , MIMO ,state-space,array of SISO,discrete time tf can be passed as a system to this function
//o/p is the cell array of the coefficients for MIMO and SISO arrays.

    [lhs,rhs]=argn(0)
    if rhs == 0 then
        error(msprintf(gettext("Not enough input arguments")))
        end
    select typeof(p)
    case "rational" then
    case "state-space" then
        p=ss2tf(p)
    else
        msprintf("\n")
        error(97,1),
    end;
    if or(size(p)<>[1 1]) then
        q=size(p);
        y=ndims(p);
        if(y==3) then
            num=cell(size(p,'r'),size(p,'c'),size(p,3));
                den=cell(size(p,'r'),size(p,'c'),size(p,3));
            for i=1:size(p,'r')
                for j=1:size(p,'c')
                    for k=1:size(p,3)
                        num(i,j,k).entries=coeff(p(i,j,k).num);
                        den(i,j,k).entries=coeff(p(i,j,k).den);
                        end
                    end
                end
            elseif(y==2) then
                //num=0;
                //den=0;
                num=cell(size(p,'r'),size(p,'c'));
                den=cell(size(p,'r'),size(p,'c'));
            for i=1:size(p,'r')
                //nn=0;
                //dd=0;
                for j=1:size(p,'c')
                    
                        num(i,j).entries= coeff(p(i,j).num);
                        den(i,j).entries= coeff(p(i,j).den);
                       // nn=[nn n];
                        //dd=[dd d];
                    end
                    //nn(1,1)=[];
                    //num=[num; nn];
                    //dd(1,1)=[];
              //      den=[den; nn];
                    end
            //num(1,:)=[];
            //den(1,:)=[];
            end
            
        else
num=coeff(p.num);//taking numerator coefficients
den=coeff(p.den);//taking denominator coefficients
//if given transfer function is continuous time then sampling time will be given as 'c'
end
if(p.dt=='d') then
    Ts=1;
else
    Ts=p.dt
end
//Ts=p.dt;
endfunction
