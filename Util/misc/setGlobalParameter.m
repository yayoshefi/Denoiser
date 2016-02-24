function setGlobalParameter()
global Parameter Analysis
%%  ############# Parameter  ##############
InitClustersNUM=350;        MaxSubSpace=0;      MinimunClusterSize=20;
SplitType='median'; % 'median' or 'totvar'
Context='comeans'; % [] or 'spectral' or 'graphcut' 'rl' or 'mrf' 'entropy'

% **************     other clustring Parameters     ****************
Parameter.MSS=MaxSubSpace;          Parameter.minclstsize=MinimunClusterSize;
Parameter.SplitType=SplitType;      Parameter.Context=Context;
Parameter.subsample=0.20;

% **************      Method Parameter     *****************
if ~isfield (Parameter,'wsize2');Parameter.wsize2=81;end  %default value
Parameter.values=struct('VarianceSplit',5*Parameter.wsize2,'Distance',InitClustersNUM,...
    'kmeans',InitClustersNUM,'Distance_Normalize',InitClustersNUM);

Parameter.values.VarianceSplit=300;
% sigma^2*wsize^2  (sigma^2=15)-lower than the noise                  cluster amount
%sigma^2*THR (THR=?5)- this is for the case of first vec energy

% **************   Spectral (KSVD) Parameters     ****************
dictsize=300;           sparsity=3;         HardThr=1;
Parameter.Spectral=struct('clustrsNUM',InitClustersNUM,'dictsize',dictsize,...
    'sparsity',sparsity,'HardThr',HardThr,'Fast',true);

% **************     Spatial clustring Parameters     ****************
Parameter.Spatil.spatialdist='decomposition';        Parameter.Spatil.lambda=0.005;
%'landmarks' / 'decomposition'  /  'simplenoramlize' /'none'
Parameter.Spatil.sigma=20;   Parameter.Spatil.NN=5;      Parameter.Spatil.CoOcThr='unused';
Parameter.location='C:/Users/Yair/Dropbox/Thesis code';
%%  ############# Analysis  ##############
Analysis=struct('Show',true,'Save',true,'figures',2,'Fast',true,'miniwindow',5,...
    'Handles',[],'K',0,'LabelsSize',...
    [Parameter.row-sqrt(Parameter.wsize2)+1,Parameter.col-sqrt(Parameter.wsize2)+1]);

end

%{
Parameter fields:

Parameter.values - is the main Paramter for each clustring Method
VarianceSplit - a paramter for max Var in a cluster
Distance & kmeans -  how much initial cluster to begin with


Parameter.Spectral - set of Parameters for spectral clustring and KSVD
'clustrsNUM'
dictsize
sparsity
HardThr



Parameter.minclstsize - cluster with less then minclstsize will be deleted
Parameter.MSS - for the case where we look for a subapce cluster
Parameter.SplitType -  in what method will the split function will work
Parameter.Context  -  which way will the spatial information will use

%}