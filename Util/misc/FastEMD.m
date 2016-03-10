function dist=FastEMD (F1,W1,F2,W2)
%% dist=FastEMD (F1,W1,F2,W2)
% FastEMD   approximated Earth Mover's Distance between two signatures
%    [dist] = FastEMD(F1,  W1, F2, W2) is the Earth Mover's Distance
%    between two signatures S1 = {F1, W1} and S2 = {F2, W2}. F1 and F2
%    consists of feature column vectors which describe S1 and S2, respectively.
%    Weights of these features are stored in W1 and W2.
%    The ground distance between two feature vectors is an euclidean 2 norm.
%    Both signatures should have the same features vectors
% 
%    for multiplie compare, use W2=row stacking of diffrent distribution
%    each row for a diffrent compare
%    the result dist in a row vector

[d,features]=size(F1);%[d,features]=size(F2)

W1=W1./sum(W1);
if isvector(W2)
    W2=W2./sum(W2);
else
    W2N=sum(W2,2);
    W2=W2./W2N(:,ones(1,size(W2,2)));
end

% COM1=sum(DiagonalMult(F1,W1','r'),2);
COM1=F1*W1';
% COM2=sum(DiagonalMult(F2,W2','r'),2);
COM2=F2*W2';

cluster1ToCOM1=W1(ones(d,1),:).*F1-COM1(:,ones(1,features));
COM1tocluster2=W2(ones(d,1),:).*F2-COM1(:,ones(1,features));

cluster2ToCOM2=W1(ones(d,1),:).*F2-COM2(:,ones(1,features));
COM2tocluster1=W2(ones(d,1),:).*F2-COM1(:,ones(1,features));


ub1= sqrt(sum( sum(cluster1ToCOM1.^2) )) + sqrt(sum( sum(COM1tocluster2.^2) ));
ub2= sqrt(sum( sum(cluster2ToCOM2.^2) )) + sqrt(sum( sum(COM2tocluster1.^2) ));

if isvector(W2)
    lb=sum(W2)*norm(COM1-COM2);
else
    lb=W2N'.*sqrt( sum(  (COM1( :,ones(size(W2,1),1) )-COM2).^2 ,1  ) );
end

%dist=(ub1+lb)/2;
dist=lb;            %by random check (curse of dimmentionalty)


end