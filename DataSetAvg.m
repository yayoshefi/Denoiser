% %%---------------- Data Collecting-------------------Git Version
clearvars 
load Data\ExpImages.mat;        load Data\rectImage.mat 
load Data\BSDS300_test.mat;     load Data\BSDS300_train.mat;

description='Avg denoising values';
%% --------------------------- PARAMETERS ------------------------------
I=BSDS300_test;         clearvars -except I description

global Parameter Analysis
Method='kmeans';        metric ='euclidean';

sigma=25;       wsize=11;       NN=9;
lambda= 0.03;   rule=3;        

MainT=tic;
Avg_Context_res=0;Avg_BM3D_res=0;Avg_ORACLE_res=0;Avg_res=0;Avg_KSVD_res=0;
Avg_Context_Sps=0;Avg_ORACLE_Sps=0;Avg_sps=0;
L=length(I);
PsnrStrct(L+1)=struct('Name',[],'Kmeans',[],'CoC',[],'ORACLE',[],'BM3D',[],'KSVD',[]);
CoOcStrct(L)=struct('Name',[],'Entropy',[],'Sparsity', [],...
        'ORACLE_Entropy',[],'ORACLE_Sparsity', [],'Context_Entropy',[],'Context_Sparsity',[]);
L=50;
for i=1:L
    Image=I(i).Image;
    name=(I(i).name);  
    if ( mod(i,5)==0 ); fprintf('  %i / %i images in %i sec\n' , i,L, round(toc(MainT)));    end;
Parameter=struct('description',description,'ImageName',name,'row',size(Image,1),'col',...
    size(Image,2),'Method',Method,'sigma',sigma,'wsize2',wsize^2,'normalize',...
    0,'metric',metric);
%% ------------------------- INITIALIZATION ---------------------------
NewPath=genpath([pwd,'/Util']);
addpath(NewPath);
setGlobalParameter();                   Analysis.Show=false;
row=Parameter.row;                      col=Parameter.col;
Parameter.spatial.NN=NN;                Parameter.spatial.UpdateRule=rule;
Analysis.DebuggerIter=1000;             Parameter.spatial.lambda=lambda;

Analysis.LabelsSize=[Parameter.row-sqrt(Parameter.wsize2)+1,Parameter.col-sqrt(Parameter.wsize2)+1];

Noise=randn(size(Image))*sigma;
Input=double(Image)+Noise;             
Data=im2col(Input,[wsize,wsize],'sliding');

X=Data;
Xmean=mean(X);
X=X-Xmean(ones(wsize^2,1),:);
Xnorm=(sum(X.^2)).^0.5;
Xn=X./(Xnorm(ones(wsize^2,1),:)+0.01);

% ORACLE
[ORACLE,ORACLECenters]=FindClusters(im2col(double(Image),[wsize,wsize],'sliding') );
[ORACLEOutput]=removenoise(double(Image),Noise,ORACLE);
ORACLEresult=psnr(ORACLEOutput,double(Image),255);
[ORACLE_Entropy, ORACLE_Sparsity]=...
    CoOc_V1 (lcm(ORACLE,Parameter.spatial.AssginType,Parameter.spatial.CoOc),false,'both',0);
% BM3D
[BM3dresult, BM3dOutput] = BM3D(im2double(Image), im2double(Image)+(Noise/255), sigma,'np',0);
%KSVD
Parameter.KSVD_params.x=Input;      Parameter.KSVD_params.sigma=sigma;
[KSVDOutput, ~] = ksvddenoise(Parameter.KSVD_params,0);
KSVDresult=psnr(KSVDOutput,double(Image),255);


if Parameter.normalize==2; Patches=Xn;
elseif Parameter.normalize==1; Patches=X; 
else Patches=Data;
end
if strcmp (Method,'gabor')  
    G=gabor ([wsize,wsize/2],[0,30,60,90,120,150]);
    Mag=imgaborfilt(Image,G);       P=(wsize-1)/2;
    Patches=(  reshape(Mag(P+1:end-P,P+1:end-P,:),[(row-2*P)*(col-2*P),length(G)] )  )';
end

%% ----------------------- Cluster -----------------------------------
[AssignVec, Centers,Energy,Basis]=...
    FindClusters(Patches,'maxsubspace',Parameter.MSS);
[Entropy,Sparsity]=...
    CoOc_V1 (lcm(AssignVec,Parameter.spatial.AssginType,Parameter.spatial.CoOc),false,'both',0);

%% ------------------ Remove Noise --------------------------------
[Output]=removenoise(double(Image),Noise,AssignVec);
result=psnr(Output,double(Image),255);
            
%% ------------------  Context  --------------------------------    
[AssignVec2]=SpatialContext (Patches,AssignVec,Centers);
[Context_Output]     =removenoise(double(Image),Noise,AssignVec2);
Context_result   = psnr(Context_Output      ,double(Image),255);

[RI,MH  ]               = RandIndex(AssignVec (:),ORACLE(:));
[Context_RI,Context_MH] = RandIndex(AssignVec2 (:),ORACLE(:));

[Context_Entropy, Context_Sparsity]=...
CoOc_V1 (lcm(AssignVec2,Parameter.spatial.AssginType,Parameter.spatial.CoOc),false,'both',0);

%% structures to save
    Psnr=struct('Name',name,'Kmeans',result,'CoC',Context_result,'ORACLE',ORACLEresult,'BM3D',BM3dresult,'KSVD',KSVDresult);
    PsnrStrct(i)=Psnr;
    CoOcs= struct('Name',name,'Entropy',Entropy,'Sparsity', Sparsity,...
        'ORACLE_Entropy',ORACLE_Entropy,'ORACLE_Sparsity', ORACLE_Sparsity,...
        'Context_Entropy',Context_Entropy,'Context_Sparsity',Context_Sparsity);
    CoOcStrct(i)=CoOcs;
%% Summary
    Avg_res=Avg_res+result;                 Avg_Context_res=Avg_Context_res+Context_result;
    Avg_BM3D_res=Avg_BM3D_res+BM3dresult;   Avg_ORACLE_res=Avg_ORACLE_res+ORACLEresult;         Avg_KSVD_res=Avg_KSVD_res+KSVDresult;
    Avg_Context_Sps=Avg_Context_Sps+Context_Sparsity;   Avg_ORACLE_Sps=Avg_ORACLE_Sps+ORACLE_Sparsity; Avg_sps=Avg_sps+Sparsity;
           
end % Image
Avg_res=Avg_res/i;  Avg_Context_res=Avg_Context_res/i;  Avg_BM3D_res=Avg_BM3D_res/i;  Avg_ORACLE_res=Avg_ORACLE_res/i;      Avg_KSVD_res=Avg_KSVD_res/i;
Avg_Context_Sps=Avg_Context_Sps/i;   Avg_ORACLE_Sps=Avg_ORACLE_Sps/i;   Avg_sps=Avg_sps/i;

PsnrStrct(i+1)=struct('Name','mean values','Kmeans',Avg_res,'CoC',Avg_Context_res,'ORACLE',Avg_ORACLE_res,'BM3D',Avg_BM3D_res,'KSVD',Avg_KSVD_res);
CoOcStrct(i+1)=struct('Name','mean','Entropy',[],'Sparsity', Avg_sps,'ORACLE_Entropy',[],...
    'ORACLE_Sparsity', Avg_ORACLE_Sps,'Context_Entropy',[],'Context_Sparsity',Avg_Context_Sps);
full_Data=struct('Psnr',PsnrStrct,'Sparsity',CoOcStrct,'Parameters',Parameter);
disp (PsnrStrct(i+1))
 %% save info
mkdir(strcat(Parameter.location,'\Results/',date));
filName=strcat('full_Data_sigma',numstr(sigma),'_clusters',num2str(Parameter.values.kmeans),'.mat');
save (strcat(Parameter.location,'\Results/',date,'/',filName), 'full_Data', '-v7.3')

