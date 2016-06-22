%%---------------- Data Collecting-------------------Git Version
clearvars ContextPsnr Psnr CentersCount Labeling new
load Database.mat;load Sport+_DB.mat; Images = createImages();
Image=barbara;      %clearvars -except Image
description='barbara Mutual';
%% --------------------------- PARAMETERS ------------------------------

global Parameter Analysis
%check struct Name
Method='kmeans';          %Distance  ,  VarianceSplit , kmeans , 'Spectral'
metric ='euclidean';        %distance function can be 'euclidean','mahalanobis'
                            %'varing_cluster_size'                            

% ---- arrays ----
sigma_array=[15,25,50,70];
wsize_array=[5,9,13,15];
normalize_array=[0];          %normalize 0-do nothing ; 1-only bias; 2- bias and gain
lambda_array=[0.001,0.01,0.05,0.1,0.3];%,0.5,0.75,0.9,0.99];%[0.01,0.05,0.1];
NN_array=[3,5,7,9,11];
CoOcThr_array=[1e-10];%[0.0005,0.0001];

Arrays=struct('Lambda',lambda_array,'NN',NN_array,'CoOcThr',CoOcThr_array,...
    'normalize',normalize_array,'wsize',wsize_array,'sigma',sigma_array);
Parameter=struct('description',description,'row',size(Image,1),'col',size(Image,2),...
    'Method',Method,'metric',metric);


%% ------------------------- INITIALIZATION ---------------------------
NewPath=genpath([pwd,'/Util']);
addpath(NewPath);
setGlobalParameter();                   Analysis.Show=true;
row=Parameter.row;                      col=Parameter.col;

Analysis.Arrays=Arrays;
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
        
        [ORACLE,ORACLECenters]=FindClusters(im2col(double(Image),[wsize,wsize],'sliding') );
        
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
                startcontext=toc;
close all
if ischar(Parameter.Context)
    [AssignVec2]=SpatialContext (Patches,AssignVec,Centers);
    contexttime=toc-startcontext;

    iter=2;
    [Context_Output]=removenoise(double(Image),Noise,AssignVec2);
    [iter_Output   ]=removenoise(double(Image),Noise,Analysis.iterations(iter).AssignVec2);
    
    result2      =psnr(Context_Output,double(Image),255);
    iter_result2=psnr(iter_Output   ,double(Image),255);

    ContextPsnr.(['Lambda',num2str(l)]).(['NN',num2str(NN)]).(['CoOcThr',num2str(thr)])...
        .(['normalize',num2str(normalize)]).(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=result2;
    
    iterPsnr.(['Lambda',num2str(l)]).(['NN',num2str(NN)]).(['CoOcThr',num2str(thr)])...
        .(['normalize',num2str(normalize)]).(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=iter_result2; 
    
    Labeling(find (wsize==wsize_array)).(['Lambda',num2str(l)]).(['NN',num2str(NN)])...
        .(['CoOcThr',num2str(thr)]).(['normalize',num2str(normalize)]).(['sigma',num2str(sigma)])=AssignVec2;

    
    fprintf(['#Centers    psnr    Context    noise    ',Parameter.Context,'-Lambda  context(m)  K2   NN   CoOc    CC_{Thr}   psnr_iter\n',...
        ' %3u       %2.3f   %2.3f     %3u            %3G           %4.3G   %u    %u    %s     %2.1G   %2.3f\n'],...
        round(size(Centers,3)), result,result2, sigma,Parameter.Spatil.lambda,...
        contexttime/60,length (unique(AssignVec2)),NN,Parameter.Spatil.CoOc,Parameter.Spatil.CoOcThr,iter_result2)

    PrintDnoise ({Output,Context_Output},{result,result2},AssignVec,AssignVec2)
    
else
    disp(['  #Centers    psnr    noise    ',Method,'    normalize',' cluster(s)'])
    disp([round(size(Centers,3)), result, round(sigma),...
        Parameter.values.(Method),Parameter.normalize ,round(itertime)])
    PrintDnoise ({Output},{result})
%     [EpsNorm1,EpsNorm2]=ShowCoOc(AssignVec);
end
% if sigma==0
%     Labeling.Oracle=AssignVec;end
% if isfield(Labeling,'Oracle')
   rand_index2=RandIndex(ORACLE(:),AssignVec2(:));
   rand_index =RandIndex(ORACLE(:),AssignVec (:));
RI.(['Lambda',num2str(l)]).(['NN',num2str(NN)]).(['CoOcThr',num2str(thr)])...
        .(['normalize',num2str(normalize)]).(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=rand_index2;
% end
    
Psnr.(['normalize',num2str(normalize)]).(['wsize',num2str(wsize)])...
    .(['sigma',num2str(sigma)])=result;
% CentersCount.(['Lambda',num2str(l)]).(['normalize',num2str(normalize)])...
%     .(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=round(size(Centers,3));
% EpsNorm.(['Lambda',num2str(l)]).(['normalize',num2str(normalize)])...
%     .(['wsize',num2str(wsize)]).(['sigma',num2str(sigma)])=[EpsNorm1,EpsNorm2];
new.(['Lambda',num2str(l)]).(['NN',num2str(NN)]).(['CoOcThr',num2str(thr)])...
        .(['normalize',num2str(normalize)]).(['wsize',num2str(wsize)])(find(sigma==sigma_array)) = ...
        struct('sigma',sigma,'PrePsnr',result,'PreRI',rand_index,'K',round(size(Centers,3)),...
        'ContextPSNR',result2,'ContextRI',rand_index2,'K2',length (unique(AssignVec2)),'iterations',Parameter.Spatil.MaxIter);

            end % Lambda
            end % NN
            end % CoOcThr
        end     %Normalize
        % ORACLE
        Labeling(find (wsize==wsize_array)).Oracle=ORACLE;
        figure;imshow(col2im(ORACLE,[wsize,wsize],[row,col]),[]);colormap jet;title (['ORACLE for w_1=',num2str(wsize)])
        Labeling(find (wsize==wsize_array)).wsize=wsize;
    end         %Wsize
end             %sigma
 %% save info
full_Data=struct('psnr',Psnr,'ContextPsnr',ContextPsnr,'IterPsnr',iterPsnr,...
    'Randindex',RI,'Parameters',Parameter,'Arrays',Arrays);
save (strcat(Parameter.location,'\Results/',date,'/','full_Data.mat'), 'full_Data', '-v7.3')
% save (strcat(Parameter.location,'\Results/',date,'/','Labeling.mat'), 'Labeling', '-v7.3')