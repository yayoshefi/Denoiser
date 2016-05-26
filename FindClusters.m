function [AssignVec, Centers,varargout]= FindClusters(Data,varargin)
%% [AssignVec, Centers,varargout]= FindClusters(Data,Parameter,varargin)
%
% Data is a n by signal_number matrix,
% Output is index matrix 1 by pnum and centers matrix 1 by 1 by K
% varargout is Basis, Energy to be used at mahalanobis ditance

%% parameter
global Parameter Analysis

pnames = {'maxsubspace'};
dflts =  {0 };
[MaxSubSpace] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});

[wsize2,pnum]=size(Data);
samp=randperm(pnum,round(Parameter.subsample*pnum));
SampData=Data(:,samp);
if Parameter.ORACLE
    if Analysis.ORACLE.level==3
        AssignVec=Analysis.ORACLE.AssignVec;    Centers=Analysis.ORACLE.Centers;
    else
    [AssignVec,~]=Dist2SubSpace(Data,Analysis.ORACLE.Centers,'dim',zeros(1,1,Parameter.values.kmeans));
    end
else
switch Parameter.Method
    case 'Distance'
        [~, Centers,Basis,Energy,E80]= DistanceClustering(SampData,MaxSubSpace);
        [AssignVec, ~]=Dist2SubSpace (Data,Centers,'basis',Basis,'dim',min(MaxSubSpace,E80));
        
    case 'VarianceSplit'
        [~, Centers,Basis,E80]= VarianceSplitClustering(SampData);
        [AssignVec, ~]=Dist2SubSpace (Data,Centers,'basis',Basis,'dim',min(Parameter.MSS,E80));
        
    case 'kmeans'
        [idx,C]=kmeans(SampData',Parameter.values.kmeans,'onlinephase','off');
        Centers=reshape(C',[size(Data,1),1,Parameter.values.kmeans]);
        [AssignVec,~]=Dist2SubSpace(Data,Centers,'dim',zeros(1,1,Parameter.values.kmeans));

    case 'Spectral'
        L=isfield(Parameter,{'clustersNUM','dictsize','sparsity','HardThr'});
        
        AssignVec=MySpectralClustering(Data,'ksvd','sparsity',...
            Parameter.Spectral.sparsity,'dictsize',Parameter.Spectral.dictsize,...
            'hardthr',Parameter.Spectral.HardThr);
        Centers=0;
        
end
end 
if nargout>2    % This can ruin the centers when using ORACLE
    [Centers,Basis,AssignVec,Energy,E80]= UpdateCenter(Data,AssignVec,false);
    varargout{1}=Energy;
    varargout{2}=Basis;
end
% Analysis.K=size(Centers,3);
Analysis.K=length(unique(AssignVec));
end