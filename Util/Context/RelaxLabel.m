function AssignVec2=RelaxLabel(Data,AssignVec,Centers)
%% AssignVec2=RelaxLabel(Data,AssignVec,Centers)
% NOT SURE WORKING...

global Parameter

K=size(Centers,3);
wsize=sqrt(Parameter.wsize2);
m=Analysis.LabelsSize(1);
n=Analysis.LabelsSize(2);

knn=0; value=true;
[E]=affinity (Data, Centers,knn,value);          %,varargin{1:6});
[E_directional,comp]=Compatability (E,AssignVec);
S=sum( E_directional.*comp(ones(K,1),:,:),3 );

p=logical(padarray(ones(m-2,n-2),[1,1]));
p=p(:)';
Enew=E(ones(K,1),p).*(1+S);
[Pr, AssignVec2]=max(Enew);

end