function Dist=HistCmp(H1,H2,method)
% Dist=HistCmp(H1,H2,method)
% method : batacheriya or chisqr
% Dist is the differnce between the to histogrmas
% if chi square: the dist is realated ti the first histogram

bins=size (H1,2);
if size (H2,2)~=bins; error ('two hist r not the same size'); end

switch method
    case 'batacheriya'
        miu1=sum(H1,2)/bins;
        miu2=sum(H2,2)/bins;
        Dist=sqrt(1- (1/sqrt(miu1*miu2*bins^2)) * sum(sqrt(H1.*H2),2) );
    case 'chisqr'
        Dist=sum( ((H1-H2).^2)./H1 ,2);
end
   
end