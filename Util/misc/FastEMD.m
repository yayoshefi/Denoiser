function dist=FastEMD (F1,W1,F2,W2)
%% dist=FastEMD (F1,W1,F2,W2)
% FastEMD   approximated Earth Mover's Distance between two signatures
%    [dist] = FastEMD(F1,  W1, F2, W2) is the Earth Mover's Distance
%    between two signatures S1 = {F1, W1} and S2 = {F2, W2}. F1 and F2
%    consists of feature column vectors which describe S1 and S2, respectively.
%    Weights of these features are stored in W1 and W2.
%    The ground distance between two feature vectors is an euclidean 2 norm.
%    Both signatures should have the same features vectors


[d,features]=size(F1);%[d,features]=size(F2)
COM1=mean(F1,2);
COM2=mean(F2,2);

cluster1ToCOM=W1(ones(d,1),:).*F1-COM1(:,ones(1,features));
COMtocluster2=W2(ones(d,1),:).*F2-COM1(:,ones(1,features));
ub= sum( sum(cluster1ToCOM.^2) ) + sum( sum(COMtocluster2.^2) );

lb=sum(W1)*norm(COM1-COM2);

dist=(sqrt(ub)+sqrt(lb))/2;


end