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
%
%    F1=[|  |  |  ...  |
%        f1 f2 f3 ...  fn
%        |  |  |  ...  |]
%
%    W1=[w1 w2 w3 ...  wn]
[d,features]=size(F1);%[d,features]=size(F2)
F1(F1==inf)=0;  F2(F2==inf)=0;      %to avoid nan option


W1=W1./sum(W1);
if isvector(W2)
    W2=W2./sum(W2);
else
    W2N=sum(W2,2);
    W2=W2./W2N(:,ones(1,size(W2,2)));
end
try
    % COM1=sum(DiagonalMult(F1,W1','r'),2);
    COM1=F1*W1';
catch ME
    msg = ['Dimension mismatch occurred: Features argument has size of: ', ...
        num2str(size(F1)),'. Weights argument has size of: ',num2str(size(W1))];
    causeException = MException('MATLAB:myCode:dimensions',msg);
    ME = addCause(ME,causeException);
    rethrow(ME);
end
% COM2=sum(DiagonalMult(F2,W2','r'),2);
COM2=F2*W2';
%{
cluster1ToCOM1=W1(ones(d,1),:).*F1-COM1(:,ones(1,features));
COM1tocluster2=W2(ones(d,1),:).*F2-COM1(:,ones(1,features));

cluster2ToCOM2=W1(ones(d,1),:).*F2-COM2(:,ones(1,features));
COM2tocluster1=W2(ones(d,1),:).*F2-COM1(:,ones(1,features));

p2p=sqrt(  sum( (W1(ones(d,1),:).*F1-W2(ones(d,1),:).*F2).^2 )  );


ub1= sqrt(sum( sum(cluster1ToCOM1.^2) )) + sqrt(sum( sum(COM1tocluster2.^2) ));
ub2= sqrt(sum( sum(cluster2ToCOM2.^2) )) + sqrt(sum( sum(COM2tocluster1.^2) ));
%}

if isvector(W2)
    lb=sum(W2)*norm(COM1-COM2);
else
    lb=W2N'.*sqrt( sum(  (COM1( :,ones(size(W2,1),1) )-COM2).^2 ,1  ) );
end

%dist=(ub1+lb)/2;
dist=lb;            %by random check (curse of dimmentionalty)
% fprintf('lower bound is: %3.2G ; Upper bpund is: %3.2G \n while the point to point dist is : %3.2G\n',lb(1),ub1(1),p2p(1));

end