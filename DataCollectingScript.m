%%---------------- Data Collecting-------------------Git Version
clearvars ContextPsnr Psnr CentersCount
load Database; %Images = createImages();
Image=barbara;
description='collecting data NN=11;Thr=0.005 NON Moving Centers ';
%% --------------------------- PARAMETERS ------------------------------

global Parameter Analysis
%check struct Name
Method='kmeans';          %Distance  ,  VarianceSplit , kmeans , 'Spectral'
metric ='euclidean';        %distance function can be 'euclidean','mahalanobis'
                            %'varing_cluster_size'                            

% ---- arrays ----
sigma_array=[10,20,30,40,60,80];
wsize_array=[5,7,9];
normalize_array=[0];          %normalize 0-do nothing ; 1-only bias; 2- bias and gain
lambda_array=[0.001,0.005,0.01];%,0.03,0.05,0.1,.3,0.5,1];
NN_array=[3,5,7,9,11];
CoOcThr_array=[0.005,0.0001];

Parameter=struct('description',description,'row',size(Image,1),'col',size(Image,2),...
    'Method',Method,'metric',metric);


%% ------------------------- INITIALIZATION ---------------------------
NewPath=genpath([pwd,'/Util']);
addpath(NewPath);
setGlobalParameter();                   Analysis.Show=true;
row=Parameter.row;                      col=Parameter.col;


for sigma=sigma_array
    Parameter.sigma=sigma;
    
    Noise=randn(size(Image))*sigma;
    Input=double(Image)+Noise;

    for wsize=wsize_array
        Parameter.wsize2=wsize^2;                                  % Parameter of varsplit is depend on wsize
        Analysis.LabelsSize=[Parameter.row-sqrt(Parameter.wsize2)+1,Parameter.col-sqrt(Parameter.wsize2)+1];
        
        Data=im2col(Input,[wsize,wsize],'sliding');
        X=Data;
        
        Xmean=mean(X);
        X=X-Xmean(ones(wsize^2,1),:);
        Xnorm=(sum(X.^2)).^0.5;
        Xn=X./(Xnorm(ones(wsize^2,1),:)+0.01);
        
        for normalize=normalize_array
        Parameter.normalize=(normalize);
        
            if Parameter.normalize==2; Patches=Xn;
            elseif Parameter.normalize==1; Patches=X; 
            else Patches=Data;
            end

%% ----------------------- Cluster -----------------------------------
            tic
            [AssignVec, Centers,Energy,Basis]=...
                FindClusters(Patches,'maxsubspace',Parameter.MSS);
            itertime=toc;
            
%% ------------------ Remove Noise --------------------------------

            [Output]=removenoise(double(Image),Noise,AssignVec);
            result=psnr(Output,double(Image),255);
            cleaningtime=toc-itertime;
            
%% ------------------  Context  --------------------------------
            for thr=1:length(CoOcThr_array);
                Parameter.Spatil.CoOcThr=CoOcThr_array(thr);
            for NN=NN_array
                 Parameter.Spatil.NN=NN;
            for l=1:length(lambda_array)
                Parameter.Spatil.lambda=lambda_array(l);
close all
if ischar(Parameter.Context)
    [AssignVec2]=SpatialContext (Patches,AssignVec,Centers);
    contexttime=toc-cleaningtime-itertime;
    [Context_Output]=removenoise(double(Image),Noise,AssignVec2);
    result2=psnr(Context_Output,double(Image),255);

    ContextPsnr.(['Lambda',num2str(l)]).(['NN',num2str(NN)]).(['CoOcThr',num2str(thr)])...
        .(['normalize',num2str(normalize)]).(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=result2;    
   
    disp(['  #Centers    psnr    Context    noise    ',Method,'    normalize',' Lambda','  context(s)'])
    disp([round(size(Centers,3)), result,result2, round(sigma),...
        Parameter.values.(Method),Parameter.normalize,Parameter.Spatil.lambda ,round(contexttime)])

    PrintDnoise ({Output,Context_Output},{result,result2},AssignVec,AssignVec2)
    
else
    disp(['  #Centers    psnr    noise    ',Method,'    normalize',' cluster(s)'])
    disp([round(size(Centers,3)), result, round(sigma),...
        Parameter.values.(Method),Parameter.normalize ,round(itertime)])
    PrintDnoise ({Output},{result})
%     [EpsNorm1,EpsNorm2]=ShowCoOc(AssignVec);
end

Psnr.(['normalize',num2str(normalize)]).(['wsize',num2str(wsize)])...
    .(['sigma',num2str(sigma)])=result;
% CentersCount.(['Lambda',num2str(l)]).(['normalize',num2str(normalize)])...
%     .(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=round(size(Centers,3));
% EpsNorm.(['Lambda',num2str(l)]).(['normalize',num2str(normalize)])...
%     .(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=[EpsNorm1,EpsNorm2];
            end % Lambda
            end % NN
            end % CoOcThr
        end
    end
end
 %% save info
full_Data=struct('psnr',Psnr,'ContextPsnr',ContextPsnr,'Parameters',Parameter);
save (strcat(Parameter.location,'\Results/',date,'/','full_Data.mat'), 'full_Data', '-v7.3')

rmpath(NewPath);