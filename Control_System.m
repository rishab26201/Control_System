VA = [-2 0 1; 1 -1 -1; -2 1 0]
VB = [1 0; 0 1; 1 0]
XA = [-1 -1 0; 1 -1 -1; 0 1 -1]
XB = [1 0; 0 0; 0 1]
T = [1 0 0; 0 0 1; 1 -1 0]

XC = [-1 0 0; 0 1 0; 0 0 1]
XD = [1 0; 0 0; 0 -1]
VC = [-1 0 0; 1 -1 0; 0 0 1]
VD = [1 0; 0 0; 0 -1]

T*XA*inv(T)
T*XB

[T1, lam1] = eigs(XA)
inv(T1)*XA*T1

[T2, lam2] = eigs(VA)

%LYAPUNOV TO BE CHECKED

syms s
I = [1 0 0; 0 1 0; 0 0 1]
XA_inv = inv(s*I - XA)
eXAt = ilaplace(XA_inv)

XB_inv = inv(s*I - T*XA*inv(T))
eXBt = ilaplace(XB_inv)


%%Lyapunov of XA
syms p1 p2 p3 p4 p5 p6
P = [p1 p2 p3; p2 p4 p5; p3 p5 p6]
mat = XA'*P + P*XA
mat = equationsToMatrix(mat)
Q = [-1; 0; 0; 0; -1; 0; 0 ;0 ;-1]
x = linsolve(mat,Q)

x = [x(1) x(2) x(3); x(2) x(4) x(5); x(3) x(5) x(6)] 


%%Lyapunov of VA
syms p1 p2 p3 p4 p5 p6
P = [p1 p2 p3; p2 p4 p5; p3 p5 p6]
mat = VA'*P + P*VA
mat = equationsToMatrix(mat)
Q = [-1; 0; 0; 0; -1; 0; 0 ;0 ;-1]
x = linsolve(mat,Q)
x = [x(1) x(2) x(3); x(2) x(4) x(5); x(3) x(5) x(6)] 

control11 = [XB XA*XB XA*XA*XB]
rank(control11)

control12 = [VB VA*VB VA*VA*VB]
rank(control12)
 
for i = 1:3
	control21 = [lam1(i)*I - XA XB]
	disp(['The rank for lamda1[' num2str(i) '] is [' num2str(rank(control21)) ']'])
end

lam1


for i = 1:3
	control22 = [lam2(i)*I - VA VB]
	disp(['The rank for lamda2[' num2str(i) '] is [' num2str(rank(control22)) ']'])
end

lam2


obs1 = [XC XA*XC XA*XA*XC]
rank(obs1)

obs2 = [VC VA*VC VA*VA*VC]
rank(obs2)