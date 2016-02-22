function [AssignVec, Centers,Basis,Energy,E80]=DistanceClustering(Data,MaxSubSpace)
%% [AssignVec, Centers,Basis,Energy,E80]=DistanceClustering(Data,Parameter,MaxSubSpace)
%

%  -----parameter-----------
global Parameter
[wsize2, pnum]=size(Data);

%% initializtion
[Centers,Basis,E100]=CreaterandCenters(wsize2,Parameter.values.Distance);
[AssignVec, tmp_dist]=Dist2SubSpace (Data,Centers,'basis',Basis,'dim',min(MaxSubSpace,E100));

switch Parameter.normalize
    case 0
        MaxRadius=norm (Data(:,1))/(wsize2^0.25);
        MaxRadius= max(tmp_dist)*0.9;   %*0.85
    case 1
        MaxRadius=norm(Data(:,1))*0.8;
        MaxRadius= max(tmp_dist);%*0.85;
    case 2
        MaxRadius=1.2*5;      %1.2*sigma*0.1 
        MaxRadius=norm(Data(:,1))*(wsize2^0.25);
end
if strcmp('mahalanobis',Parameter.metric)
    MaxRadius=0.45;  %0.42
end
disp (['initializing Cluster maximum radius to: ',num2str(MaxRadius)])
ChangedAssign=pnum;
LastAssignment=AssignVec;
KChange=zeros(2,5);
clearvars tmp_dist;

iter=0;
%% main loop
while ChangedAssign>pnum*0.01
    iter=iter+1;
    if ~mod(iter,5); DeBugger();end

    K=size(Centers,3);
    [Centers,Basis,~,E,E80]= UpdateCenter(Data,AssignVec);
        
    [AssignVec, Distances]=Dist2SubSpace (Data,Centers,'basis',Basis,'dim',min(MaxSubSpace,E80),'s',E);
    ChangedAssign=sum(AssignVec~=LastAssignment);

    KChange(1,:)=[KChange(1,2:end),K-size(Centers,3)];
    
    LargeClusters=unique( AssignVec(Distances>MaxRadius) );

    [Centers, Basis,AssignVec,E80]=...
        SplitClusters(Data,Centers,Basis,LargeClusters,AssignVec,E80);   
        
    LastAssignment=AssignVec;
    ChangeRadius()
     KChange(2,:)=[KChange(2,2:end),size(Centers,3)-K];
    
    
end
[Centers,Basis,AssignVec,Energy,E80]= UpdateCenter(Data,AssignVec);

%% --------------    Util Functions   -----------------------
    function [Centers,Basis,E100]=CreaterandCenters(dim1,dim3)

    Centers=Data(:,randperm(pnum,dim3));
    Centers=reshape(Centers,[dim1,1,dim3]);

    Basis=repmat(eye(dim1),1,1,dim3);

    E100=dim1*ones(1,1,dim3);
    end

%{
    function deleteNAN()
        Centers(isnan(Centers))=[]; Centers=reshape(Centers,wsize2,1,[]);
        Basis(isnan(Basis))=[];     Basis=reshape(Basis,wsize2,wsize2,[]);
        E(isnan(E))=[];             E=reshape(E,wsize2,1,[]);
        E80(isnan(E80))=[];         E80=reshape(E80,1,1,[]);
    end
%}

    function  DeBugger()
        disp (['iter: ',num2str(iter)]); 
        disp ([num2str(ChangedAssign), ' points changed assignment'])
%         disp ([num2str(CenterMovment), ' Center Movment '])
%         disp ([num2str(maxMovment), ' Maximum Center Movment '])
        if sum(find(KChange))
            disp ([num2str(KChange(1,:)), ' deleted Centers last 5 iterations'])
            disp ([num2str(KChange(2,:)), ' Change  Centers last 5 iterations'])
        end
        disp (['total ',num2str(size(Centers,3)), ' Centers in this iteation'])
        disp (['max distance  ', num2str(max(Distances))])
%         n=mod(iter/5,4)-1;
%         if ~n; figure; end
%         subplot(2,2,n+2);title (['tmp lables, iterate: ',num2str(iter)])
%         imagesc(col2im(AssignVec,[sqrt(Parameter.wsize2),sqrt(Parameter.wsize2)],[Parameter.row,Parameter.col],'sliding'));
        
    end

    function ChangeRadius()
        if size(Centers,3)<Parameter.values.Distance*0.8        %if # Clusters is too small
            str='Decreased max Radius: ';
            MaxRadius=MaxRadius*0.9;
        elseif size(Centers,3)>Parameter.values.Distance*1.3
            str='Enlarged max Radius ';
            if ~isempty(LargeClusters)
                str=strcat(str,' (',num2str(numel(LargeClusters)),' Large clusters): ');
                MaxRadius=MaxRadius*1.25;
            end
        else    return;
            
        end
        disp (['iteration #',num2str(iter),' ',str,num2str(round(MaxRadius))])
    end

end
% 0) set random centers and Basis, and find assignment
% 1) update centers (mean of cluster) and Basis (svd to cluster), delete
%       empty clusters if necessry
% 2) split the very big clusters
% 3) find closest SubSpace {center and subspace} to each data point to decide clusters
% 4) if  data points still change assignment cluters, go to -> 1