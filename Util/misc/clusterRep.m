function [C,H]=clusterRep(Data,AssignVec,Hist)
%% [C,H]=clusterRep(Data,AssignVec,,Hist)
% this function find cluster Represeantative in both visual space and Hist
% space. the function is used in mutual dist context method.
% 
% C- the visual Space rep. is merely the mean of all patchs in cluster
% C(:,:,j) is the Represeantative(Center) of cluster j
%
% H [K by K]- the Hist Space Hist, is one of the histograms in the cluster
% which has the smallest sum dist to all other hist in the cluster.
% H(j,:) is the Represeantative of cluster j.
% 

[Centers,~,AssignVec,~,~]=UpdateCenter(Data,AssignVec,false);
K=size(Centers,3); pnum=size(Data,2);

H=zeros(K,K);
tic
for k=1:K
    cluster=Hist(:,AssignVec==k);
    clstsize=size(cluster,2);
    RepD=inf; Rep=zeros(1,K);
    for c_i=1:clstsize
        candidate=cluster(:,c_i);
        tempD=0;
        for m=1:clstsize
            tempD=tempD+emd(squeeze(Centers)',squeeze(Centers)',candidate,cluster(:,m),@gdf);
        end
        if tempD<RepD
            RepD=tempD;  Rep=candidate;   end
    end
    H(k)=Rep;
    if ~mod(k,10); disp (strcat('calculated: ',num2str(k), ' clusters Hist Rep. in :',num2str (toc), ' sec'));end
end
C=Centers; 
end