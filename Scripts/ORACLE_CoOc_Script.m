% %%---------------- CoOc ORACLE -------------------Git Version
p=pwd;      NewPath=genpath(p(1:end-8));    addpath(NewPath);
clearvars 
load Data\ExpImages.mat;        load Data\rectImage.mat 
load Data\BSDS300_test.mat;     load Data\BSDS68.mat;

description='ORACLE the CoOc matrix';
%% --------------------------- PARAMETERS ------------------------------
I=Irect;         clearvars -except I description

global Parameter Analysis
Method='kmeans';        metric ='euclidean';

sigma=50;      wsize=11;       CoOcType='CC';
rule=3;        beta=0.005;
lambda_array= 3E-4;%[10.^[-6:-2] 5*10.^[-5:-4] ];   %[0:0.05:0.2 0.4:0.2:0.8];


for i=12
    Image=I(i).Image;
    name=(I(i).name);   disp(name) 
    
Parameter=struct('description',description,'ImageName',name,'row',size(Image,1),'col',...
size(Image,2),'Method',Method,'sigma',sigma,'wsize2',wsize^2,'normalize',...
0,'metric',metric);

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

%% Adding Noise 
Patches=Data;
K=size(ORACLECenters,3);              M=length(Data);
OEACLEDistances=zeros (K,size(Data,2));
for k=1:K
    pointdist=sum(Data.^2)-2*ORACLECenters(:,:,k)'*Data+ORACLECenters(:,:,k)'*ORACLECenters(:,:,k);
    OEACLEDistances(k,:)=abs(pointdist).^0.5;
end
ind=sub2ind([K,M],ORACLEVec,1:M);
AvgORACLEDist=sum(OEACLEDistances(ind) )/M;

%% Clustering with Noise 
[AssignVec, Centers,~,~] = FindClusters(Patches,'maxsubspace',Parameter.MSS);
[SimpleOutput]=removenoise(double(Image),Noise,AssignVec);
Simpleresult=psnr(SimpleOutput,double(Image),255);
PsnrStrct(1)= struct('Params',[name, ': K-means'],'Psnr',Simpleresult,'improve',ORACLEresult-Simpleresult,...
    'Score',0,'iterations',[],'dist',[],'l_1',[],'Labels',[],'K',[]);

%% Learned CoOc
 Analysis.ORACLE=struct('Centers',ORACLECenters,'CoOc',ORACLECoOc,'AssignVec',ORACLEVec,'level',2);
 Parameter.ORACLE=true;
 
[AssignVec, Centers,~,~] = FindClusters(Patches,'maxsubspace',Parameter.MSS);
CoOc = lcm(AssignVec,Parameter.CoOc.AssginType,Parameter.CoOc.Type);
[Entropy,l_0,l_1] = CoOc_V1 (CoOc,false,'both',Parameter.CoOc.epsilon);
ind=sub2ind([K,M],AssignVec,1:M);
AvgDist=sum(OEACLEDistances(ind) )/M;

% D noise
[Output]=removenoise(double(Image),Noise,AssignVec);
result=psnr(Output,double(Image),255);

fprintf ('ORACLE denoising, Psnr=%2.3f , margin (%1.4f)\n',ORACLEresult,ORACLEresult-result)
PsnrStrct(2)= struct('Params',['K-means','; level=1 '],'Psnr',result,'improve',ORACLEresult-result,...
    'Score',AvgDist + beta*l_1,'iterations',[],'dist',AvgDist,'l_1',l_1,...
    'Labels',reshape( uint16(AssignVec) ,Analysis.LabelsSize),'K',length(unique(AssignVec)));

for NN=3
    Parameter.spatial.NN=NN;
for l=1:length(lambda_array)
    Parameter.spatial.lambda=lambda_array(l);
            
[AssignVec2]        = SpatialContext (Patches,AssignVec,Centers);
[Context_Output]    = removenoise(double(Image),Noise,AssignVec2);
Context_result      = psnr(Context_Output      ,double(Image),255);

ContextCoOc = lcm(AssignVec2,Parameter.CoOc.AssginType,Parameter.CoOc.Type);
[Context_Entropy, Context_l_0,Context_l_1] = CoOc_V1 (ContextCoOc ,false,'both',Parameter.CoOc.epsilon);

% fprintf ('Using: lambda=%G , NN=%i\n \t score=%3.1f Psnr=%2.3f (%1.4f) \n',...
%     lambda_array(l),NN,score,Context_result,Context_result-result)

ind=sub2ind([K,M],AssignVec2',1:M);
AvgContext_Dist=sum(OEACLEDistances(ind) )/M;
score = AvgContext_Dist + beta * Context_l_1;

if isfield(Analysis,'iterations')
        iter=Analysis.iterations;
else    iter=[];    end
PsnrStrct(end+1)=struct('Params',['lambda=',num2str(lambda_array(l),'%1.3G'),'; NN=',num2str(NN)],...
    'Psnr',Context_result,'improve',Context_result-result,'Score',score,'iterations',iter,...
    'dist',AvgContext_Dist,'l_1',Context_l_1,'Labels',reshape( uint16(AssignVec2) ,Analysis.LabelsSize),...
    'K',length(unique(AssignVec2)));

end %lambda
end %NN
end %image
PsnrStrct(end+1)=struct('Params','ORACLE','Psnr',ORACLEresult,'improve',0,...
    'Score',AvgORACLEDist + beta*ORACLE_l_1,'iterations',[],'dist',AvgORACLEDist,'l_1',ORACLE_l_1,...
    'Labels',reshape( uint16(ORACLEVec) ,Analysis.LabelsSize),'K',length(unique(ORACLEVec)));

%% sumary
full_Data=struct('Psnr',PsnrStrct,'Parameters',Parameter);
T = struct2table(rmfield(PsnrStrct,{'iterations','Labels'}));
disp(T)

%  %% save info
% mkdir(strcat(Parameter.location,'\Results/',date));
% filName=strcat('full_Data_sigma',num2str(sigma),'_clusters',num2str(Parameter.values.kmeans),'.mat');
% save (strcat(Parameter.location,'\Results/',date,'/',filName), 'full_Data', '-v7.3')



