//Function:- DC gain///////////////////
//cont time transfer function SISO
s=poly(0,'s');
q=syslin('c',(s+8)/(s^2+8*s+4));
z=dcgain(q);
//MIMO//
s=poly(0,'s');
sysa=syslin('c',(s+2)/(s^2+2*s+3))
sys1=syslin('c',(s+8)/(s^8+3*s+4));
sysm=[sysa sys1];
z1=dcgain(sysm);
//SISO array//
aa=pid(rand(2,3,3),3,4,5);
z2=dcgain(aa);
//State-space//
bb=ssrand(2,2,3);
z3=dcgain(bb);
