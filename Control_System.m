VA = [-2 0 1; 1 -1 -1; -2 1 0]
VB = [1 0; 0 1; 1 0]
XA = [-1 -1 0; 1 -1 -1; 0 1 -1]
XB = [1 0; 0 0; 0 1]
T = [1 0 0; 0 0 1; 1 -1 0]

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

