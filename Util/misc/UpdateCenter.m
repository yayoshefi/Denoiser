function [Centers,Basis,AssignVec,Energy,E80]=UpdateCenter(Data,AssignVec,varargin)
%% [Centers,Basis,AssignVec,ENERGY,E80]=UpdateCenter(Data,AssignVec,varargin)
% Updates the cluster Center and Basis according to current Assignment
%
% varargin is bollean for:  true= assign '0' to small clusters
%                           false= do nothing to small clusters
%
% in the case a cluster is too small it deletes that cluster
% Energy is(wsize2-by 1-by K) a column vectors with variance of each basis vectors
% for the normalize case, the bias Vec Variance is set to 1.
% E80 is (1-by-1 -by -K) what siae of the Basis contains 80 precent Energy

% Parameters
global Parameter

K=max(AssignVec);
[wsize2, pnum]=size(Data);
asgn.type='()';

Centers=zeros(wsize2,1,K);
tmpE=zeros(1,wsize2);
E80=zeros(1,1,K);

Basis=zeros(wsize2,wsize2,K);       %full size Basis
Energy=zeros(wsize2,1,K);

if nargin==2 varargin{1}=true; end

for k=K:-1:1
    clstsize=sum(AssignVec==k);
    if clstsize<Parameter.minclstsize && varargin{1}
        Centers(:,:,k)=[];   Basis(:,:,k)=[];
        Energy(:,:,k)=[];          E80(:,:,k)=[];
        AssignVec(AssignVec==k)=0;
        AssignVec=AssignVec-(AssignVec>k);  %updates the allready done centers
        continue;       %skipes that Data, it will have no center
    end            
    cluster=Data(:,AssignVec==k);
    Centers(:,1,k)=mean(cluster,2);
    if Parameter.normalize==2
        Centersnorm=(sum(Centers(:,1,k).^2)).^0.5;
        Centers(:,1,k)=Centers(:,1,k)./Centersnorm(ones(wsize2,1),:,:);
    end


    cluster_0=cluster-Centers(:,ones(1,clstsize),k);
    [U,S]=svd(cluster_0,'econ');

    S=sum(S.^2,2);
    asgn.subs={1:size(S,1)};
    E=subsasgn(tmpE,asgn,S);
    Energy(:,:,k)=E';

    E=cumsum(E);
    Eper=E/E(end);
    if clstsize<2
        E80(k)=0;
    else
        E80(k)=find(Eper>0.8,1);
    end

    asgn.subs={1:size(U,1),1:size(U,2),k};
    Basis=subsasgn(Basis,asgn,U);


end
end