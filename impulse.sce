//execute function pole,pid first
//Function:-impulse
s=poly(0,'s');
sys=syslin('c',(s+3)/(s^3+4*s+2));
impulse(sys)

sys1=ssrand(2,3,4);
impulse(sys1)

impulse(sys,sys1)

impulse(sys,'--r',sys1,'gx')

aa=pid(rand(2,3,4),2,3,4);
impulse(aa,%T)

[y t x ysd]=impulse(sys1);
[y t]=impulse(sys);
