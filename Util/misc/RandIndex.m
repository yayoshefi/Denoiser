function RI=RandIndex(L1,L2)
%RANDINDEX - calculates Rand index to compare two partitions

K=max(L1);
M=max(L2);
N=length (L1);

A=0;        % # pairs that are in same      cluster in L1 & L2
B=0;        % # pairs that are in diffrent  cluster in L1 % L2
for k=1:K
    a=zeros(1,M);
    for m=1:M
        agreements= (L1==k) & (L2==m);
        a(m)=sum(agreements);      % points that r in same clusters in 2 domains
        if a(m)>1
            A=A+nchoosek(a(m),2);
        end
        L1not=(L1~=k);        L2not=(L2~=m);
        OtherLabel=L1not & L2not;   % points that r not in either clusters
        B=B+a(m)*sum(OtherLabel);
    end
    %duplicates??
end
B=B/2;
RI=(A+B)/(nchoosek(N,2));
end
