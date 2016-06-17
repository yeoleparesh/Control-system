



//Function:- bandwidth////////////////////////
//state-space SISO with default db drop=-3db
a1=ssrand(1,1,2);
y=bandwidth(a1);
///with given dbdrop
y1=bandwidth(a1,-2);
//tf model//
s=poly(0,'s');
sys=syslin('c',(s+2)/(s^4+3*s+12));
y2=bandwidth(sys);
//discrete time
sys=syslin(0.2,(s+2)/(s^4+3*s+12));
y3=bandwidth(sys);
//
//SISO array --default dbdrop ("when passed through pid value of bandwidth will be nan or inf")
aa=pid(rand(2,2),3,4,5);
y4=bandwidth(aa,[],1);
//SISO array with given dbdrop
y5=bandwidth(aa,-4,1);





