function setGlobalParameter()
%%  USING ORACLE
%   ORACLE=struct(Centers,CoOc,AssignVec,level)
%   ORACLE Level:
%       1- is only using centers
%       2- is using centers an CC
%       3- is for comparing only and gives the true 0 noise labeling
%
%   To create the ORACLE run main with:
%   no context no noise and normalize=-5 
%   this will create an ORACLE of the Assignment and Centers.
%   Build the CoOc matrix using ShowCoOc(AssignVec)
%
%
%   To Use The defined ORACLE:
%       SET  Parameter.ORACLE=true
%   it will create the field Analysis.ORACLE and update it there exist a
%   variable named ORACLE
%
global Parameter Analysis ORACLE
%%  ############# Parameter  ##############
InitClustersNUM=250;        MaxSubSpace=0;      MinimunClusterSize=20;
Debug=false;                USEORACLE=false;    ORACLE_level=1;
SplitType='median'; % 'median' or 'totvar'
Context='comeans'; % [] or 'spectral' or 'graphcut' 'rl' or 'mrf' 'entropy'
                   % 'mutualdist', 'comeans' 
% **************     other clustring Parameters     ****************
Parameter.MSS=MaxSubSpace;          Parameter.minclstsize=MinimunClusterSize;
Parameter.SplitType=SplitType;      Parameter.Context=Context;
Parameter.subsample=0.20;           Parameter.ORACLE=USEORACLE;

% **************      Method Parameter     *****************
if ~isfield (Parameter,'wsize2');Parameter.wsize2=81;end  %default value
Parameter.values=struct('VarianceSplit',5*Parameter.wsize2,'Distance',InitClustersNUM,...
    'kmeans',InitClustersNUM,'Distance_Normalize',InitClustersNUM);
Parameter.values.VarianceSplit=300;

% **************   Spectral (KSVD) Parameters     ****************
dictsize=300;           sparsity=3;         HardThr=1;
Parameter.Spectral=struct('clustrsNUM',InitClustersNUM,'dictsize',dictsize,...
    'sparsity',sparsity,'HardThr',HardThr,'Fast',true);

% **************     Context clustring Parameters     ****************
Parameter.Spatil.spatialdist='decomposition';        Parameter.Spatil.lambda=0.75;
%'landmarks' / 'decomposition'  /  'simplenoramlize' /'none'
Parameter.Spatil.sigma=20;   Parameter.Spatil.NN=5;      Parameter.Spatil.CoOcThr='unused';
Parameter.location='C:/Users/Yair/Dropbox/Thesis code';
%%  ############# Analysis  ##############
if Parameter.ORACLE; if ~exist('ORACLE','var');ORACLE=Analysis.ORACLE;end;end %save last ORACLE Or update
Analysis=struct('Show',true,'Save',true,'figures',1,'Fast',true,'Handles',[]...
    ,'K',0,'DebuggerMode',Debug,'ColorMap','lines','LabelsSize',...
    [Parameter.row-sqrt(Parameter.wsize2)+1,Parameter.col-sqrt(Parameter.wsize2)+1]);
if Parameter.ORACLE;
    Analysis.ORACLE=ORACLE; Analysis.ORACLE.level=ORACLE_level;end
end
%%  USING ORACLE
%   first use no context no noise and normalize=-5 
%   -> get the ORACLE for AssignVec and for Centers. Build CoOc
%   create structure: ORACLE=struct(Centers,CoOc,AssignVec)
%   run the algorithm with Parameter.ORACLE=ture and it will use the built
%   ORACLE

%%
% Parameter.values.VarianceSplit=300;
% % sigma^2*wsize^2  (sigma^2=15)-lower than the noise                  cluster amount
% %sigma^2*THR (THR=?5)- this is for the case of first vec energy

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