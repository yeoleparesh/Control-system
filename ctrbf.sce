//Function ctrbf

//Example(1):-
A=[1 1;4 -2];B=[1 -1;1 -1];C=[1 0;0 1];
[a b c k]=ctrbf(A,B,C)
//with tolerance

[at bt ct kt]=ctrbf(A,B,C,10)

[at bt ct kt t]=ctrbf(A,B,C,10)

//Example(2):-
L=[1 2;3 4];M=[2 3;4 5];N=[3 4;4 5];O=[3 4;5 6];
[p q r]=ctrbf(L,M,N)
 
//references:-
// https://in.mathworks.com/help/control/ref/ctrbf.html
