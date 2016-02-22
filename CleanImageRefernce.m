%%---------------- Clean Image Refernce-------------------Octover Version
% clear all;%clc;close all;
load Database
Image=barbara;
%%--------------------------- PARAMETERS ------------------------------

Method='Distance';   %check struct Name        %Distance  ,  VarianceSplit , kmeans 'Distance_Normalize'
Parameter.normalize=2;       %normalize 0-do nothing ; 1-only bias; 2- bias and gain
Parameter.metric ='euclidean';  %distance function can be 'euclidean','mahalanobis'


clustersNUM=300;        Parameter.MaxSubSpace=0;
dictsize=300;           sparsity=1;     HardThr=1;
varparmeter=1;

%% ------------------------- INITIALIZATION ---------------------------
[row,col]=size(Image);

for wsize=9
    Data=im2col(double(Image),[wsize,wsize],'sliding');
    X=Data;

    Xmean=mean(X);
    X=X-Xmean(ones(wsize^2,1),:);
    Xnorm=(sum(X.^2)).^0.5;
    Xn=X./(Xnorm(ones(wsize^2,1),:)+0.01);
    
    if Parameter.normalize==2; Patches=Xn;
    elseif Parameter.normalize==1; Patches=X; 
    else Patches=Data;
    end
    Parameter.VarianceSplit=170; Parameter.kmeans=clustersNUM;Parameter.Distance=clustersNUM;
    Parameter.Spectral=[clustersNUM,dictsize,sparsity,HardThr];
    Parameter.Distance_Normalize=clustersNUM;
    

    %% ----------------------- Cluster -----------------------------------
    
    tic
    [AssignVec, Centers]=...
        FindClusters(Patches,Method,Parameter.(Method),...
        'metric',Parameter.metric,'maxsubspace',Parameter.MaxSubSpace);
    itertime=toc;

    %% ------------------ Remove Noise --------------------------------
    for sigma=20:20:100

        Input=double(Image)+randn(size(Image))*sigma;

        [Output]=removenoise(Input,wsize,AssignVec,sigma);
        result=psnr(Output,double(Image),255);

        Psnr.(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=result;
%         Pics.(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=Output;
        CentersCount.(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=round(size(Centers,3));
        
        disp(['  #Centers    psnr    noise    ',Method,'    clusteringtime'])
        disp([round(size(Centers,3)), result, round(sigma), Parameter.(Method) ,round(itertime)])
    end
end

%% ------------------- Save info ------------------------------------
Clean_DistaneREGULAR_Nor_2_subspace_0=struct('psnr',Psnr,'centers',CentersCount,'parametes',Parameter);
save Clean_DistaneREGULAR_Nor_2_subspace_0.mat Clean_DistaneREGULAR_Nor_2_subspace_0
