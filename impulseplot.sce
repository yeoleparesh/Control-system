//Function:-impulseplot

s=poly(0,'s');
sys=syslin('c',(s+3)/(s^3+4*s+2));
h=impulse(sys)
set(h.children,"foreground",13);
set(h.children,"polyline_style",2);  
sys1=ssrand(2,3,4);
impulse(sys1)
impulse(sys,sys1)
impulse(sys,'--r',sys1,'gx')
aa=pid(rand(2,3,4),2,3,4);
impulse(aa,%T)
