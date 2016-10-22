// Function sigma
// Example1 
s=%s;
sys=syslin('c',1/(1+s));
sigma(sys)
sigma(sys,sys^2)
sigma(sys,20:40)
[V,w]=sigma(sys)
