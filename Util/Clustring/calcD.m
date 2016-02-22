function d = calcD (E)
%% d = calcD (E) 
%
% INPUT:
% E is [n by K] matrix, such that, Sim= E*E' 
% we calaulte it using deconposition of
% columns for E and 'rows' of E', using the fact
% that Sim is  symetric D can also be sum of columns
%
% OUTPUT:
% d is [1Xn] vector  which can populate 
% the D diagonal matrix
% calculate  Diaginal Volume matrix
% without using the Full similarity matrix


[n,K]=size (E);      % n is the num of signals

%% Memory conservation
% C=zeros(1,K);
% d=zeros(1,n);
% for chunk=1:10000:n
%     C=C+( sum(abs( E(:,chunk:min(chunk+9999,n)) ),2) )';
% end
% 
% for chunk=1:10000:n
%     d(chunk:min(chunk+9999,n))=C*abs(E(:,chunk:min(chunk+9999,n)));
% end 
%% -------------------------------------------------------------

if min (E) < 0
    E=abs(E);end

C=sum(E);
d=C*E';

% D=sparse (1:n,1:n,d,n,n);

end