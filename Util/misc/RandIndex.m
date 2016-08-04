function [RI,MH]=RandIndex(L1,L2)
%RANDINDEX - calculates Rand index to compare two partitions
L1=L1(:); L2=L2(:);

P=max(L1);
Q=max(L2);
N=length (L1);

A=0;        % # pairs that are in same      cluster in L1 & L2
B=0;        % # pairs that are in diffrent  cluster in L1 % L2
for p=1:P
    a=zeros(1,Q);
    for q=1:Q
        agreements= (L1==p) & (L2==q);
        a(q)=sum(agreements);      % points that r in same clusters in 2 domains
        if a(q)>1
            A=A+nchoosek(a(q),2);
        end
        L1not=(L1~=p);        L2not=(L2~=q);
        OtherLabel=L1not & L2not;   % points that r not in either clusters
        B=B+a(q)*sum(OtherLabel);
        
        M(p,q)=a(q);          % Confiusion matrix - comparing clustering - jan\2007
    end
    %duplicates??
end
B=B/2;
RI=(A+B)/(nchoosek(N,2));

MH=mean(max(M,[],2));     % L2 is the ORACLE
end


