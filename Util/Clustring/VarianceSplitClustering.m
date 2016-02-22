function [AssignVec, Centers,Basis,E80]=VarianceSplitClustering(Data,MaxSubSpace)
%% [AssignVec, Centers,Basis,E80]=VarianceSplitClustering(Data,Parameter,metric,MaxSubSpace)
% Var split is a Tree like algorithm that keeps spliting the clusters until
% each cluster has a small enough variance

%% initialization
global Parameter
[wsize2, pnum]=size(Data);
AssignVec=ones(1,pnum);
LargeClusters=1;
[Centers,Basis,AssignVec,E,E80]=UpdateCenter(Data,AssignVec,Parameter);

iter=0;
while ~isempty(LargeClusters)
    
    iter=iter+1;
%     if ~mod(iter,10)
%         disp (['iter: ',num2str(iter)])
%         disp ([num2str(length(LargeClusters)), ' cluster r still too spread'])
%         disp ([num2str(size(Centers,3)), ' Centers in this iteation'])
%         disp (['max cluster Energy  ', num2str(max( sum(E(:,end,:)) ))])
%         disp (['max distance  ', num2str(max(Distances))])
%         disp (['min distance  ', num2str(min(Distances))])
%         
%     end
    

    [Centers, Basis,AssignVec,E80]=...
        SplitClusters(Data,Centers,Basis,LargeClusters,AssignVec,E80);
    
%     [Centers,Basis,AssignVec,E,E80]=UpdateCenter(Data,AssignVec);
    
    [AssignVec, Distances]=Dist2SubSpace (Data,Centers,'basis',Basis,Parameter,'dim',min(MaxSubSpace,E80),'s',E);
    
    [Centers,Basis,AssignVec,E,E80]=UpdateCenter(Data,AssignVec,Parameter);
        
    LargeClusters=find( sum(E(:,end,:))>Parameter.values.VarianceSplit );
    
    
    
end

end


% 0) set zero center and set all asiignVec to him 
% 1) split large clusters (median of cluster_0- using the center)
% 2) update2center and find Centers
% 3)compute cluster properties: largest var (first find assignment Dist2subspace)
% large cluster are when the total energy is bigger than parmeter (2 option is only first main vec)
% 4) if large cluster still exisit go to-> 1