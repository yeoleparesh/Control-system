//Function:-impulseplot

s=poly(0,'s');
sys=syslin('c',(s+3)/(s^3+4*s+2));
h=impulseplot(sys)
set(h.children,"foreground",13);
 
sys1=ssrand(2,3,4);
impulseplot(sys1)

impulseplot(sys,sys1)

impulseplot(sys,'--r',sys1,'gx')

aa=pid(rand(2,3,4),2,3,4);
impulseplot(aa,%T)
