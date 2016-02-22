function D=PseudoDiagInv (d)
%% D=PseudoDiagInv (d)
% d is a diagonal represantation, can be a mtrix or an vec of the main diag
% D=d^-1 for diagonal mtrix
% D is the same size of d

L_m=ismatrix(d);  

if L_m
    d=diag(d); end

d(abs(d)<10^6*eps)=1;
D=1./d;

if L_m
    D=diag(D); end

end