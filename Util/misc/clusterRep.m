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

%% Visual space
[Centers,~,AssignVec,~,~]=UpdateCenter(Data,AssignVec,false);
K=size(Centers,3); pnum=size(Data,2);

%% Histogrma Space
HistSpace='Lb_EMD';
switch HistSpace
    case 'CC'
        H=ShowCoOc(AssignVec,'CoOc');close (gcf);

    case 'Lb_EMD'
        H=zeros(K,K);
        strt=toc;
        for k=1:K
            cluster=Hist(AssignVec==k,:);
            clstsize=size(cluster,1);
            RepD=inf; Rep=zeros(1,K);
            for c_i=1:clstsize
                candidate=cluster(c_i,:);
                tempD=FastEMD(squeeze(Centers),candidate,squeeze(Centers),cluster);
                tempD=sum(tempD);
        %         tempD=0;        
        %         for m=1:clstsize
        %             tempD=tempD+emd(squeeze(Centers)',squeeze(Centers)',candidate',cluster(:,m)',@gdf);
        %         end
                if tempD<RepD
                    RepD=tempD;  Rep=candidate;   end
            end
            H(k,:)=Rep;
            samp=20;
            if ~mod(k,samp); fprintf ('calculated total: %u Histograms clusters Rep;   Last %u clusters in :% 4.2g sec\n',k,samp,toc-strt);strt=toc;end
        end
        fprintf('found representatives\n');
end
C=Centers; 
end