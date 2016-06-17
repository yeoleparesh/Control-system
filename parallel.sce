//Function:-parallel
//Example taken from --https://books.google.co.in/books?id=RGPNBQAAQBAJ&pg=PA480&lpg=PA480&dq=connect+parallel+function+matlab&source=bl&ots=Ll9TA0DJKP&sig=PnZZVZFQ7TDttAe7MdcgFmjVx2M&hl=en&sa=X&ved=0ahUKEwjStZ_moKzNAhUGsJQKHTV3BJIQ6AEIVjAJ#v=onepage&q=connect%20parallel%20function%20matlab&f=false
s=%s;sys1=syslin('c',1/(500*s^2));sys2=syslin('c',(s+1)/(s+2));
sys=parallel(sys1,sys2);

//
 
a1=[1 2;3 4];b1=[2 3;4 5];c1=[3 4;5 6];d1=[4 5;6 7];a2=[1 4;3 4];b2=[2 6;4 5];c2=[7 4;5 6];d2=[4 9;6 7];
aa=syslin('c',a1,b1,c1,d1);bb=syslin('c',a2,b2,c2,d2);
sys=parallel(aa,bb,1,2,2,1);

[a b c d]=parallel(a1,b1,c1,d1,a2,b2,c2,d2,1,2,2,1);
