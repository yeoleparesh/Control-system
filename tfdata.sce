////Function:-- tfata ////////////
///cont. transfer function model//
s=poly(0,'s');
sysa=syslin('c',(s+2)/(s^2+2*s+3));
[num,den,Ts]=tfdata(sysa);
///discrete transfer function model////
z=poly(0,'z');
sys=syslin(0.02,(s+2)/(s^2+2*s+3));
[num1,den1,Ts1]=tfdata(sys);
///MIMO//
sys1=syslin('c',(s+8)/(s^8+3*s+4));
sysm=[sysa sys1];
[num2,den2,Ts2]=tfdata(sysm);
////SISO array///////
aa=pid(rand(2,3,3),3,4,5);
[num3,den3,Ts3]=tfdata(aa);
////state-space///
bb=ssrand(2,2,3);
[num4,den4,Ts4]=tfdata(bb);
