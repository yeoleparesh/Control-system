// Function sigmaplot
// Example1 
s=%s;
sys=syslin('c',1/(1+s));
sigmaplot(sys)
sigmaplot(sys,sys^2)
sigmaplot(sys,20:40)
