% %%---------------- CoOc ORACLE -------------------Git Version
clearvars 
load Data\ExpImages.mat;        load Data\rectImage.mat 
load Data\BSDS300_test.mat;     load Data\BSDS68.mat;

description='ORACLE the CoOc matrix';
%% --------------------------- PARAMETERS ------------------------------
I=I;         clearvars -except I description

global Parameter Analysis
Method='kmeans';        metric ='euclidean';

sigma=50;      wsize=11;       CoOcType='MI';
rule=3;        beta=0.1;
lambda_array= [0.001];


for i=7:7
    Image=I(i).Image;
    name=(I(i).name);   disp(name) 
    
Parameter=struct('description',description,'ImageName',name,'row',size(Image,1),'col',...
size(Image,2),'Method',Method,'sigma',sigma,'wsize2',wsize^2,'normalize',...
0,'metric',metric);
NewPath=genpath([pwd,'/Util']);
addpath(NewPath);
setGlobalParameter();

Parameter.spatial.CoOc=CoOcType;            Parameter.spatial.UpdateRule=rule;
Analysis.DebuggerIter=5;
    
Noise=randn(size(Image))*sigma;
Input=double(Image)+Noise;             
Data=im2col(Input,[wsize,wsize],'sliding');

X=Data;
Xmean=mean(X);
X=X-Xmean(ones(wsize^2,1),:);
Xnorm=(sum(X.^2)).^0.5;
Xn=X./(Xnorm(ones(wsize^2,1),:)+0.01);

%% ORACLE
[ORACLEVec,ORACLECenters]=FindClusters(im2col(double(Image),[wsize,wsize],'sliding') );
ORACLECoOc=lcm(ORACLEVec,Parameter.spatial.AssginType,Parameter.spatial.CoOc);
[ORACLE_Entropy, ORACLE_Sparsity] = CoOc_V1 (ORACLECoOc,false,'both',0);
[ORACLEOutput]=removenoise(double(Image),Noise,ORACLEVec);
ORACLEresult=psnr(ORACLEOutput,double(Image),255);

 Analysis.ORACLE=struct('Centers',ORACLECenters,'CoOc',ORACLECoOc,'AssignVec',ORACLEVec,'level',2);
 Parameter.ORACLE=true;

%% Clustering with Noise    
Patches=Data;
[AssignVec, Centers,Energy,Basis]=...
    FindClusters(Patches,'maxsubspace',Parameter.MSS);
CoOc = lcm(AssignVec,Parameter.spatial.AssginType,Parameter.spatial.CoOc);
[Entropy,Sparsity]=...
    CoOc_V1 (CoOc,false,'both',0);
% D noise
[Output]=removenoise(double(Image),Noise,AssignVec);
result=psnr(Output,double(Image),255);
fprintf ('ORACLE denoising, Psnr=%2.3f , margin (%1.4f)\n',ORACLEresult,ORACLEresult-result)

for NN=3:2:9
    Parameter.spatial.NN=NN;
for l=1:length(lambda_array)
    Parameter.spatial.lambda=lambda_array(l);
            
[AssignVec2]=SpatialContext (Patches,AssignVec,Centers);
[Context_Output]     =removenoise(double(Image),Noise,AssignVec2);
Context_result   = psnr(Context_Output      ,double(Image),255);

ContextCoOc = lcm(AssignVec2,Parameter.spatial.AssginType,Parameter.spatial.CoOc);
[Context_Entropy, Context_Sparsity] = CoOc_V1 (ContextCoOc ,false,'both',0);

score= Analysis.iterations(end).AvgDist+beta*Analysis.iterations(end).epsNorm;

fprintf ('Using: lambda=%G , NN=%i\n \t score=%3.1f Psnr=%2.3f (%1.4f) \n',...
    lambda_array(l),NN,score,Context_result,Context_result-result)

end %lambda
end %NN
end %image

%% sumary
% full_Data=struct('Psnr',PsnrStrct,'Sparsity',CoOcStrct,'Parameters',Parameter);
% disp (PsnrStrct(i+1))
%  %% save info
% mkdir(strcat(Parameter.location,'\Results/',date));
% filName=strcat('full_Data_sigma',num2str(sigma),'_clusters',num2str(Parameter.values.kmeans),'.mat');
% save (strcat(Parameter.location,'\Results/',date,'/',filName), 'full_Data', '-v7.3')



