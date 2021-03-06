function [AssignVec, Distances]=Dist2SubSpace (Data,Centers,varargin)
%% [AssignVec, Distances]=Dist2SubSpace (Data,Centers,varargin)
%
% Input: Data is [dim by pnum], each column is a data point in dim space
% Centers is [dim by 1 by K], each center is in dim space
%
% Output
%  AssignVec is [knn by pnum], the Label set for each Data point.
%  When knn>1 the Assignvec is sorted by Distances to each Label
%
%  Distances is [knn by pnum], the distance between each Data point
%  and the corospondibg Centers in AssignVec
%  Distances(1,:)- is the distance from each data to the closest Center
%
% [AssignVec, Distances]=Dist2SubSpace (...,'basis',B).         default{0}
% [AssignVec, Distances]=Dist2SubSpace (...,'dim',E80).         default{0}
% [AssignVec, Distances]=Dist2SubSpace (...,'s',E).             default{1}
% E is [dim by 1 by K] the autovariance of each Center component
% [AssignVec, Distances]=Dist2SubSpace (...,'knn',knn).         default{1}
% knn- K Nearest Neighbour
%
% the function uses Parameters set in globar variable 'Parameter':
%
%   metric: 'euclidean'{default} or 'mahalanobis'
%       euclidean norm2 Distace is computed: ||(p-c)||2-||B'(p-c)||2
%       when E80=0, we compute the Distance ||(p-c)||2, with no Basis B
%       mahalanobis : divide each compomnent by its standad devation

% Parameters
global Parameter
K=size(Centers,3);
[wsize2, pnum]=size(Data);
pnames = {'basis'   'dim'        's'                'knn'};
dflts =  {false,zeros(1,1,K), ones(wsize2,1,K),   1};       %changed Basis from zeros(wsize2,wsize2,K) because of memory issues
[Basis,Dim,S,knn] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});

knn=min(K,knn);
AssignVec=zeros(knn,pnum);
Distances=inf*ones(knn,pnum);
E=inf*ones(K,pnum);

S=PseudoDiagInv(S);



h=waitbar(0,'Computing data distances..','CreateCancelBtn','setappdata(gcbf,''cancel'',1)');
setappdata(h,'cancel',0);
switch Parameter.metric
    case 'euclidean'
        for k=1:K
            
            if getappdata(h,'cancel')
            break
            end
            waitbar(k/K);
            if Parameter.wsize2==1;
                pointdist=Data.^2-2*Centers(:,:,k)'*Data+Centers(:,:,k)'*Centers(:,:,k);
            else
                pointdist=sum(Data.^2)-2*Centers(:,:,k)'*Data+Centers(:,:,k)'*Centers(:,:,k);
            end
            if Centers(1,:,k)==inf; pointdist=inf; end %fix empty clusters
        %     pointdist=sum( (Data-Centers(:,ones(1,pnum),k)).^2 );
           if Dim(k)==0
               projectiondist=0;
           else
               projectiondist=sum( (Basis(:,1:Dim(k),k)'*Data).^2 );
           end
           tmpdist=(abs(pointdist-projectiondist)).^0.5;
           if knn==1
               track=tmpdist<Distances(end,:);
               Distances(end,track)=tmpdist(track);
               AssignVec(end,track)=k;
           else
               E(k,:)=tmpdist;
           end
        end
           if knn>1
               waitbar(k/(K+10),h,'Sorting smallest distances ...');
               [Distances,I]=sort(E);
               Distances=Distances(1:knn,:);
               
               %sort the matrix AssignVec by distances
               [col_sub,AssignMat]=meshgrid(1:pnum,1:K);
               sort_indx=sub2ind([K,pnum],I(:),col_sub(:));
               AssignVec=reshape(AssignMat(sort_indx),K,[]);
               AssignVec=AssignVec(1:knn,:);

           end

    case 'mahalanobis'
        for k=1:K
            
            if getappdata(h,'cancel')
            break
            end
            waitbar(k/K);
            
%             var=sum(S(:,:,k).^2,2);
%             pointdist=sum(Data.^2)-2*Centers(:,:,k)'*Data+Centers(:,:,k)'*Centers(:,:,k);
            pointdist=sum( ((Basis(:,:,k)'*Data-Basis(:,:,k)'*Centers(:,ones(1,pnum),k)).^2)...
                .*S(1:end,ones(1,pnum),k) );
            tmpdist=(abs(pointdist)).^0.5;
            if knn==1
                track=tmpdist<Distances(end,:);
                Distances(end,track)=tmpdist(track);
                AssignVec(end,track)=k;
            else
                E(k,:)=tmpdist;
            end
        end
           if knn>1
                waitbar(k/K+10,h,'Sorting smallest distances ...');
                [Distances,I]=sort(E);
                Distances=Distances(1:knn,:);
               
                [col_sub,AssignMat]=meshgrid(1:pnum,1:K);
                sort_indx=sub2ind([K,pnum],I(:),col_sub(:));
                AssignVec=reshape(AssignMat(sort_indx),K,[]);
                AssignVec=AssignVec(1:knn,:);
           end
           
% ##########################################################           
% ###########  Spatial Distance metric #####################
% ##########################################################
    case 'regularized euclidean'
        lambda=0.5;
        miniwindow=5;
        for k=1:K
            pointdist=sum(Data.^2)-2*Centers(:,:,k)'*Data+Centers(:,:,k)'*Centers(:,:,k);
        %     pointdist=sum( (Data-Centers(:,ones(1,pnum),k)).^2 );
           if Dim(k)==0
               projectiondist=0;
           else
               projectiondist=sum( (Basis(:,1:Dim(k),k)'*Data).^2 );
           end
           tmpdist=(abs(pointdist-projectiondist)).^0.5;
%            NNdist=lambda*(miniwindow^2-)
           
           if knn==1
               track=tmpdist<Distances(end,:);
               Distances(end,track)=tmpdist(track);
               AssignVec(end,track)=k;
           else
               E(k,:)=tmpdist;
           end
        end
           if knn>1
               [Distances,I]=sort(E);
               Distances=Distances(1:knn,:);
               
                [col_sub,AssignMat]=meshgrid(1:pnum,1:K);
                sort_indx=sub2ind([K,pnum],I(:),col_sub(:));
                AssignVec=reshape(AssignMat(sort_indx),K,[]);
                AssignVec=AssignVec(1:knn,:);

           end
    case 'varing_cluster_size'
        lambda=0.9;             %1 unit will be between (lambda , 1/lambda)
        Parameter.metric='euclidean';
        [~, Alpha]=Dist2SubSpace (squeeze(Centers),ones(wsize2,1));
        Alpha1=min(Alpha);      Alpha2=max(Alpha);
        m=( (lambda^-1) - lambda )/(Alpha2-Alpha1);
        a=lambda-m*Alpha1;
        S=a+m*Alpha;
        % s=m*x + (Alpah2*lamda-Alpah1*lamda^-1)/(Alpha2-Alpha1)
        for k=1:K
            pointdist=sum( ((Data-Centers(:,ones(1,pnum),k)).^2) * S(k) );
            if Dim(k)==0
               projectiondist=0;
            else
               projectiondist=sum( (Basis(:,1:Dim(k),k)'*Data).^2 );
            end
           tmpdist=(abs(pointdist-projectiondist)).^0.5;
          
           if knn==1
               track=tmpdist<Distances(end,:);
               Distances(end,track)=tmpdist(track);
               AssignVec(end,track)=k;
           else
               E(k,:)=tmpdist;
           end
        end
           if knn>1
               [Distances,I]=sort(E);
               Distances=Distances(1:knn,:);
               
                [col_sub,AssignMat]=meshgrid(1:pnum,1:K);
                sort_indx=sub2ind([K,pnum],I(:),col_sub(:));
                AssignVec=reshape(AssignMat(sort_indx),K,[]);
                AssignVec=AssignVec(1:knn,:);

           end
           Parameter.metric='varing_cluster_size';
           
end
delete(h);

end
%(||p||^2-||B'p||^2)^1/2 

