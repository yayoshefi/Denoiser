function setGlobalParameter()
%%  USING ORACLE
%   global ORACLE  - so it can be read at this function
%   ORACLE=struct(Centers,CoOc,AssignVec,level)
%   ORACLE Level:
%       0- is only for comparesens
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
InitClustersNUM=50;            Context='comeans';      UpdateRule=3; 
Debug=false;    DebugIter=10;   USEORACLE=false;      ORACLE_level=2;
SplitType='median'; % 'median' or 'totvar'
%Context options:[]\'spectral'\'graphcut'\'rl'\'mrf' 'entropy' 'mutualdist' 'comeans'
MaxSubSpace=0;              MinimunClusterSize=20;
% **************     other clustring Parameters     ****************
Parameter.MSS=MaxSubSpace;          Parameter.minclstsize=MinimunClusterSize;
Parameter.SplitType=SplitType;      Parameter.Context=Context;
Parameter.subsample=0.20;           Parameter.ORACLE=USEORACLE;

% **************     Context clustring Parameters     ****************
Parameter.spatial.spatialdist='decomposition';          Parameter.spatial.lambda=0.08;
%'landmarks' / 'decomposition'  /  'simplenoramlize' /'none'
Parameter.spatial.sigma=20;   Parameter.spatial.NN=3;   Parameter.spatial.CoOcThr='unused';
Parameter.spatial.UpdateRule=UpdateRule;                Parameter.spatial.CoOc='MI';
                                                        Parameter.spatial.AssginType='hard';

Parameter.location='C:/Users/Yair/Dropbox/Thesis code';
% **************      Method Parameter     *****************
if ~isfield (Parameter,'wsize2');Parameter.wsize2=81;end  %default value
Parameter.values=struct('VarianceSplit',5*Parameter.wsize2,'Distance',InitClustersNUM,...
    'kmeans',InitClustersNUM,'Distance_Normalize',InitClustersNUM);
Parameter.values.VarianceSplit=300;
% **************   Spectral (KSVD) Parameters     ****************
dictsize=300;           sparsity=3;         HardThr=1;
Parameter.Spectral=struct('clustrsNUM',InitClustersNUM,'dictsize',dictsize,...
    'sparsity',sparsity,'HardThr',HardThr,'Fast',true);

%%  ############# Analysis  ##############
if Parameter.ORACLE; if isstruct('ORACLE');ORACLE=Analysis.ORACLE;end;end %save last ORACLE Or update
Analysis=struct('Show',true,'Save',true,'figures',1,'Fast',true,'Handles',[]...
    ,'K',0,'DebuggerMode',Debug,'DebuggerIter',DebugIter,'ColorMap','jet','LabelsSize',...
    [Parameter.row-sqrt(Parameter.wsize2)+1,Parameter.col-sqrt(Parameter.wsize2)+1]);
if Parameter.ORACLE;
    fprintf ('You choose to use ORACLE with level %u!\n',ORACLE_level)
    Analysis.ORACLE=ORACLE; Analysis.ORACLE.level=ORACLE_level;end
%% ############ Warnings ##################
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','stats:kmeans:FailedToConverge');
end

%%
% Parameter.values.VarianceSplit=300;
% % sigma^2*wsize^2  (sigma^2=15)-lower than the noise                  cluster amount
% %sigma^2*THR (THR=?5)- this is for the case of first vec energy
