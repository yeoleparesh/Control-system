////Function:-- pole ///
//transfer function model////////
s=poly(0,'s');
sysa=syslin('c',(s+2)/(s^2+2*s+3));
p=pole(sysa)
z=poly(0,'z');
sys=syslin(0.02,(s+2)/(s^2+2*s+3));
p1=pole(sys)
//MIMO system or cell array 
sys1=syslin('c',(s+8)/(s^8+3*s+4));
sysm=[sysa sys1];
r=pole(sysm);
