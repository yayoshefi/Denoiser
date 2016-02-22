function [AssignVec, Centers,Basis,Energy,E80]=DistanceClustering_Normalize(Data,Parameter,metric,MaxSubSpace)
%% [AssignVec, Centers,Basis,Energy,E80]=DistanceClustering_Normalize(Data,Parameter,metric,MaxSubSpace)
% this function gets Un-Normalized data, and computes distance using 1Vec,
% such that the 'bias' normilazation is done inside this function.
%parameter
normalize=0;
MaxRadius=450; %1.2*5;      %1.2*sigma*0.1 
if strcmp('mahalanobis',metric)
    MaxRadius=0.45;  %0.42
end

%% initializtion
[wsize2, pnum]=size(Data);
[Centers,Basis,E100]=CreaterandCenters(wsize2,Parameter);

LastCenters=Centers;        %1)movment
CenterMovment=Parameter;

LastAssignment=zeros(1,pnum);
% [Basis,~]=svd(Data,'econ');
[AssignVec, Distances]=Dist2SubSpace (Data,Centers,Basis,'metric','euclidean',...
    'dim',min(MaxSubSpace,E100),'normalizedata',normalize);
MaxRadius=max(Distances)/1.4;  %%  FIX (vary a lot between large and small noises)

ChangedAssign=pnum;

iter=0;
%% main loop
while ChangedAssign>pnum*0.005                         %CenterMovment>1
    iter=iter+1;
    if ~mod(iter,5)
        disp (['iter: ',num2str(iter)])
        disp ([num2str(ChangedAssign), ' points changed assignment'])
%         disp ([num2str(CenterMovment), ' Center Movment '])
%         disp ([num2str(maxMovment), ' Maximum Center Movment '])
        disp ([num2str(size(Centers,3)), ' Centers in this iteation'])
        disp (['max distance  ', num2str(max(Distances))])
    end

    
    [Centers,Basis,AssignVec,E,E80]=...
        UpdateCenter(Data,AssignVec,'normalizedata',normalize);
%     CenterMovment=mean( (sum( (Centers-LastCenters).^2 )).^0.5 );
%     maxMovment=max( sum( (Centers-LastCenters).^2 ).^0.5 );

        
    [AssignVec, Distances]=Dist2SubSpace (Data,Centers,Basis,'metric',metric,...
        'dim',min(MaxSubSpace,E80),'s',E,'normalizedata',normalize);
    ChangedAssign=sum(AssignVec~=LastAssignment);
    
    LargeClusters=unique( AssignVec(Distances>MaxRadius) );

    [Centers, Basis,AssignVec,E80]=SplitClusters(Data,Centers,Basis,...
        LargeClusters,AssignVec,E80,'normalizedata',normalize,'maxsubspace',MaxSubSpace);
    
    if size(Centers,3)<Parameter*0.8        %if # Clusters is too small
        MaxRadius=MaxRadius*0.9;         % keep spliting
    elseif size(Centers,3)>Parameter*1.3
        if ~isempty(LargeClusters)
        MaxRadius=MaxRadius*1.25;end
    end
    
    LastCenters=Centers;
    LastAssignment=AssignVec;
    
end
[Centers,Basis,AssignVec,Energy,E80]=...
    UpdateCenter(Data,AssignVec,'normalizedata',normalize);
end

function [Centers,Basis,E100]=CreaterandCenters(dim1,dim3)

Y=255*rand(dim1,1,dim3);

% Ymean=mean(Y);
% Y=Y-Ymean(ones(dim1,1),:,:);
% Ynorm=(sum(Y.^2)).^0.5;
% Centers=Y./Ynorm(ones(dim1,1),:,:);
Centers=Y;

Basis=repmat([eye(dim1),(1/(dim1^0.5))*ones(dim1,1)],1,1,dim3);

E100=dim1*ones(1,1,dim3);
end

% 0) set random centers and Basis, and find assignment
% 1) update centers (mean of cluster) and Basis (svd to cluster), delete
%       empty clusters if necessry
% 2) split the very big clusters
% 3) find closest SubSpace {center and subspace} to each data point to decide clusters
% 4) if  data points still change assignment cluters, go to -> 1