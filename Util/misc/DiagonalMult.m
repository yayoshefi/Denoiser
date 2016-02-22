function B=DiagonalMult(M,D,D_direction)
%% function B=DiagonalMult(M,D,direction)
% D is a Col vector [n X 1]
% B is the result of multiple M by the diagonal D
%
% B=M*D (direction='r')
% column j in M is multiplied by D(j)
%
% or
%
% B=D*M (direction='l')
% row i in M is multiplied by D(i)
[m,n]=size(M);
% B=sparse([],[],[],m,n,nz);
if strcmp(D_direction,'r')&& (m>100000)
    band=800;
else
    band=5000;
end
switch D_direction
    case 'l'                    %B=D*M:  every row is multiple of the diagonal
        for bp=1:band:m
            chunk=bp:min(bp+band-1,m);
            B(chunk,:)=M(chunk,:).*D(chunk,ones(1,n));
        end
    case 'r'                    %B=M*D   every col is multiple of the diagonal
        for bp=1:band:n
            chunk=bp:min(bp+band-1,n);
            B(:,chunk)=M(:,chunk).*D(chunk,ones(m,1))';
        end
end

end