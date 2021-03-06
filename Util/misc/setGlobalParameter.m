function setGlobalParameter()
%%  USING ORACLE
%   global ORACLE  - so it can be read at this function
%   ORACLE=struct(Centers,CoOc,AssignVec,level)
%   ORACLE Level:
%       0- is only for comparesens  = Do nothing
%       1- Uses centers
%       2- Uses both Co-Occurrence matrix and the centers 
%       3- ORACLE labeling. uses the pre-defiend 0 noise clusterting
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
InitClustersNUM=300;            Context='comeans';      UpdateRule=3;   Type='CC'; 
Debug=false;    DebugIter=5;   USEORACLE=false;      ORACLE_level=0;
SplitType='median'; % 'median' or 'totvar'
%Context options:[]\'spectral'\'graphcut'\'rl'\'mrf' 'entropy' 'mutualdist' 'comeans'
MaxSubSpace=0;              MinimunClusterSize=20;
% **************     other clustring Parameters     ****************
Parameter.MSS=MaxSubSpace;          Parameter.minclstsize=MinimunClusterSize;
Parameter.SplitType=SplitType;      Parameter.Context=Context;
Parameter.subsample=0.20;           Parameter.ORACLE=USEORACLE;

% **************     Context clustring Parameters     ****************
Parameter.spatial.spatialdist='decomposition';          Parameter.spatial.lambda=0.0000001;
%'landmarks' / 'decomposition'  /  'simplenoramlize' /'none'
Parameter.spatial.sigma=20;             Parameter.spatial.NN=9;   
Parameter.spatial.UpdateRule=UpdateRule;                

Parameter.CoOc=struct('Type',Type,'AssginType','hard','ShrinkType','matrix',...
    'ShrinkPer',0.0,'Thr',0,'epsilon',setEpsilon(Type),'logarithmic',[]);

loc=cd;     Parameter.location=strrep([loc( 1:strfind(loc,'\Doc') ),'Dropbox\Thesis code'],'\','/');
% **************      Method Parameter     *****************
if ~isfield (Parameter,'wsize2');Parameter.wsize2=81;end  %default value
Parameter.values=struct('VarianceSplit',5*Parameter.wsize2,'Distance',InitClustersNUM,...
    'kmeans',InitClustersNUM,'Distance_Normalize',InitClustersNUM,'gabor',InitClustersNUM);
Parameter.values.VarianceSplit=300;
% **************   De-Nosier Parameters     ****************
Parameter.DeNoise=struct('Robust',false,'soft',false);                       %Check
% **************   Spectral (KSVD) Parameters     ****************
dictsize=300;           sparsity=3;         HardThr=1;
Parameter.Spectral=struct('clustrsNUM',InitClustersNUM,'dictsize',dictsize,...
    'sparsity',sparsity,'HardThr',HardThr,'Fast',true);
Parameter.KSVD_params=struct('x',[],'blocksize',8,'dictsize',256,'sigma',[],'maxval',255,...
                             'trainnum',40000,'iternum',20,'memusage','normal');

%%  ############# Analysis  ##############
if Parameter.ORACLE; if isstruct('ORACLE');ORACLE=Analysis.ORACLE;end;end %save last ORACLE Or update
Analysis=struct('Show',false,'Save',true,'figures',1,'Fast',true,'Handles',[]...
    ,'K',0,'DebuggerMode',Debug,'DebuggerIter',DebugIter,'ColorMap','jet','LabelsSize',...
    [Parameter.row-sqrt(Parameter.wsize2)+1,Parameter.col-sqrt(Parameter.wsize2)+1],'samp',[]);
if Parameter.ORACLE;
    fprintf ('You choose to use ORACLE with level %u!\n',ORACLE_level)
    Analysis.ORACLE=ORACLE; Analysis.ORACLE.level=ORACLE_level;end
%% ############ Warnings ##################
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','stats:kmeans:FailedToConverge');
end
%% De-Bugging
%   DebuggerMode - sets to save iterations structure
%   DebuggerIter - when smaller than MaxIter: prints progress on command window
%   Show         - sets to print figures of P_i and P_Ni


%%
% Parameter.values.VarianceSplit=300;
% % sigma^2*wsize^2  (sigma^2=15)-lower than the noise                  cluster amount
% %sigma^2*THR (THR=?5)- this is for the case of first vec energy
