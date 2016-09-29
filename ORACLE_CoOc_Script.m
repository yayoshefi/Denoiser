% %%---------------- CoOc ORACLE -------------------Git Version
clearvars 
load Data\ExpImages.mat;        load Data\rectImage.mat 
load Data\BSDS300_test.mat;     load Data\BSDS68.mat;

description='ORACLE the CoOc matrix';
%% --------------------------- PARAMETERS ------------------------------
I=Irect;         clearvars -except I description

global Parameter Analysis
Method='kmeans';        metric ='euclidean';

sigma=50;      wsize=11;       CoOcType='CC';
rule=3;        beta=0.1;
lambda_array= [10.^[-6:-3] 5*10.^[-6:-3] ];   %[0:0.05:0.2 0.4:0.2:0.8];


for i=12
    Image=I(i).Image;
    name=(I(i).name);   disp(name) 
    
Parameter=struct('description',description,'ImageName',name,'row',size(Image,1),'col',...
size(Image,2),'Method',Method,'sigma',sigma,'wsize2',wsize^2,'normalize',...
0,'metric',metric);
NewPath=genpath([pwd,'/Util']);
addpath(NewPath);
setGlobalParameter();

Parameter.CoOc.Type=CoOcType;            Parameter.spatial.UpdateRule=rule;
Analysis.DebuggerIter=50;    Analysis.Show=false;        Analysis.DebuggerMode=true;
setEpsilon ();

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
ORACLECoOc=lcm(ORACLEVec,Parameter.CoOc.AssginType,Parameter.CoOc.Type);
[ORACLE_Entropy,ORACLE_l_0,ORACLE_l_1] = CoOc_V1 (ORACLECoOc,false,'both',Parameter.CoOc.epsilon);
[ORACLEOutput]=removenoise(double(Image),Noise,ORACLEVec);
ORACLEresult=psnr(ORACLEOutput,double(Image),255);

%% Clustering with Noise 
Patches=Data;

[AssignVec, Centers,~,~] = FindClusters(Patches,'maxsubspace',Parameter.MSS);
[SimpleOutput]=removenoise(double(Image),Noise,AssignVec);
Simpleresult=psnr(SimpleOutput,double(Image),255);
PsnrStrct(1)= struct('Params',[name, ': K-means'],'Psnr',Simpleresult,'improve',ORACLEresult-Simpleresult,'Score',0,'iterations',[]);

%% Learned CoOc
 Analysis.ORACLE=struct('Centers',ORACLECenters,'CoOc',ORACLECoOc,'AssignVec',ORACLEVec,'level',2);
 Parameter.ORACLE=true;
 
[AssignVec, Centers,~,~] = FindClusters(Patches,'maxsubspace',Parameter.MSS);
CoOc = lcm(AssignVec,Parameter.CoOc.AssginType,Parameter.CoOc.Type);
[Entropy,Sparsity,l_1] = CoOc_V1 (CoOc,false,'both',Parameter.CoOc.epsilon);
% D noise
[Output]=removenoise(double(Image),Noise,AssignVec);
result=psnr(Output,double(Image),255);

fprintf ('ORACLE denoising, Psnr=%2.3f , margin (%1.4f)\n',ORACLEresult,ORACLEresult-result)
PsnrStrct(2)= struct('Params',['K-means','; level=1 '],'Psnr',result,'improve',ORACLEresult-result,'Score',0,'iterations',[]);

for NN=9
    Parameter.spatial.NN=NN;
for l=1:length(lambda_array)
    Parameter.spatial.lambda=lambda_array(l);
            
[AssignVec2]=SpatialContext (Patches,AssignVec,Centers);
[Context_Output]     =removenoise(double(Image),Noise,AssignVec2);
Context_result   = psnr(Context_Output      ,double(Image),255);

ContextCoOc = lcm(AssignVec2,Parameter.CoOc.AssginType,Parameter.CoOc.Type);
[Context_Entropy, Context_Sparsity,Context_l_1] = CoOc_V1 (ContextCoOc ,false,'both',0);

score= Analysis.iterations(end).AvgDist+beta*Analysis.iterations(end).l_0;

% fprintf ('Using: lambda=%G , NN=%i\n \t score=%3.1f Psnr=%2.3f (%1.4f) \n',...
%     lambda_array(l),NN,score,Context_result,Context_result-result)

if isfield(Analysis,'iterations')
        iter=Analysis.iterations;
else    iter=[];    end
PsnrStrct(end+1)=struct('Params',['lambda=',num2str(lambda_array(l),'%1.3G'),'; NN=',num2str(NN)],...
    'Psnr',Context_result,'improve',Context_result-result,'Score',score,'iterations',iter);

end %lambda
end %NN
end %image

%% sumary
full_Data=struct('Psnr',PsnrStrct,'Parameters',Parameter);
T = struct2table(rmfield(PsnrStrct,'iterations'));
disp(T)

%  %% save info
% mkdir(strcat(Parameter.location,'\Results/',date));
% filName=strcat('full_Data_sigma',num2str(sigma),'_clusters',num2str(Parameter.values.kmeans),'.mat');
% save (strcat(Parameter.location,'\Results/',date,'/',filName), 'full_Data', '-v7.3')



