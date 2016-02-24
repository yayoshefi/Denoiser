%%---------------- Data Collecting-------------------Octover Version
clearvars ContextPsnr Psnr CentersCount
load Database; %Images = createImages();
Image=barbara;
description='test for rmf';
%% --------------------------- PARAMETERS ------------------------------

global Parameter Analysis
%check struct Name
Method='kmeans';          %Distance  ,  VarianceSplit , kmeans , 'Spectral'
metric ='euclidean';        %distance function can be 'euclidean','mahalanobis'
                            %'varing_cluster_size'                            

% ---- arrays ----
sigma_array=[40];
wsize_array=9;
normalize_array=[0,1];          %normalize 0-do nothing ; 1-only bias; 2- bias and gain
lambda_array=[0,0.1,1];

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
        Parameter.wsize=wsize;                                  % Parameter of varsplit is depend on wsize
        
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
            for l=1:length(lambda_array)
                Parameter.Spatil.lambda=lambda_array(l);

if ischar(Parameter.Context)
    [AssignVec2]=SpatialContext (Patches,AssignVec,Centers);
    [Context_Output]=removenoise(double(Image),Noise,AssignVec2);
    result2=psnr(Context_Output,double(Image),255);

    ContextPsnr.(['Lambda',num2str(l)]).(['normalize',num2str(normalize)])...
        .(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=result2;    
   
    disp(['  #Centers    psnr    Context    noise    ',Method,'    normalize',' cluster (s)'])
    disp([round(size(Centers,3)), result,result2, round(sigma),...
        Parameter.values.(Method),Parameter.normalize ,round(itertime)])

    PrintDnoise ({Output,Context_Output},{result,result2},AssignVec,AssignVec2)
else
     PrintDnoise ({Output},{result},AssignVec)
end

Psnr.(['Lambda',num2str(l)]).(['normalize',num2str(normalize)])...
    .(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=result;
CentersCount.(['Lambda',num2str(l)]).(['normalize',num2str(normalize)])...
    .(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=round(size(Centers,3));
close all
            end
        end
    end
end
 %% save info
Distance_Normalzie=struct('psnr',Psnr,'centers',CentersCount,'Parameters',Parameter);
save (strcat('Results/',date,'/','Distance_Normalzie.mat'), 'Distance_Normalzie', '-v7.3')
%     Parameter.location= strcat('Results/',date,'/',str);

rmpath(NewPath);