%%---------------- Main Script: De noising-------------------GitHub Version
% clearvars
load Database;Images = createImages();
Image=barbara;
description='EMD VS FastEMD';
%%--------------------------- PARAMETERS ------------------------------
global Parameter Analysis

Method='kmeans';        %Distance  ,  VarianceSplit , kmeans , 'Spectral'
sigma=40;
wsize=7;
normalize=0;            %normalize 0-do nothing ; 1-only bias; 2- bias and gain (-5)- Oracle
metric ='euclidean';    %distance function can be 'euclidean','mahalanobis'
                        %'varing_cluster_size'

Parameter=struct('description',description,'row',size(Image,1),'col',...
    size(Image,2),'Method',Method,'sigma',sigma,'wsize2',wsize^2,'normalize',...
    normalize,'metric',metric);


%% ------------------------- INITIALIZATION ---------------------------
NewPath=genpath([pwd,'/Util']);
addpath(NewPath);
setGlobalParameter();

row=Parameter.row;      col=Parameter.col;
Noise=randn(size(Image))*sigma;
Input=double(Image)+Noise;
Data=im2col(double(Input),[wsize,wsize],'sliding');
X=Data;

% X=g*p+b       NORMALIZATION
Xmean=mean(X);
X=X-Xmean(ones(wsize^2,1),:);               % X- bias normalize
Xnorm=(sum(X.^2)).^0.5;
Xn=X./(Xnorm(ones(wsize^2,1),:)+0.001);     % Xn- gain normalize

% if Analysis.UseOracle;          Patches=im2col(double(Image),[wsize,wsize],'sliding'); %OracleMode
if Parameter.normalize==2;      Patches=Xn;
elseif Parameter.normalize==1;  Patches=X; 
elseif Parameter.normalize==0;  Patches=Data;
elseif Parameter.normalize==-5; Patches=im2col(double(Image),[wsize,wsize],'sliding'); %OracleMode
end
clearvars Xmean Xnorm row col normalize description metric 
clearvars lena boat house barbara mond mondrian be
%% ----------------------- Cluster -----------------------------------
tic
[AssignVec, Centers,Energy,Basis]= FindClusters(Patches,'maxsubspace',Parameter.MSS);
itertime=toc;

%% ------------------ Remove Noise --------------------------------
[Output]=removenoise(double(Image),Noise,AssignVec);
result=psnr(Output,double(Image),255);
cleaningtime=toc-itertime;

%% -------------------------------------------------------------
%% ------------------  Context  --------------------------------
%% -------------------------------------------------------------
if ischar(Parameter.Context)
%     if Analysis.UseOracle;Patches=Data;end  %use noisy Image with no normalize
    [AssignVec2]=SpatialContext (Patches,AssignVec,Centers);

    [Context_Output]=removenoise(double(Image),Noise,AssignVec2);

    result2=psnr(Context_Output,double(Image),255);
    Out={Output,Context_Output};
    res={result,result2};
else  AssignVec2=[];
end
contexttime=toc-cleaningtime-itertime;

%% ----------------- Visualization ----------------------------
if Analysis.Show
    ShowClusters(Image,Noise,AssignVec,Centers,Energy);
    if ischar(Parameter.Context)
        ShowClusters(Image,Noise,AssignVec2,Centers,Energy);  %Centers & Energy are not updated
        PrintDnoise (Out,res);
    else PrintDnoise ({Output},{result});   
    end
else
    if ischar(Parameter.Context); PrintDnoise (Out,res);
    else PrintDnoise ({Output},{result});
    end
end

visualtime=toc-contexttime-cleaningtime-itertime;


%% ----------------------  Summary  ------------------
fprintf(['#Centers    psnr    noise    ',Method,'    normalize',' cluster(s)','  Clean(s)','  visual(s)\n',...
    ' %3u       %2.3f   %3u      %u            %u      %3G     %3G     %3G\n'],...
    size(Centers,3), result, sigma,Parameter.values.(Method),Parameter.normalize ,itertime, cleaningtime,visualtime);

% disp(['  #Centers    psnr    noise    ',Method,'    normalize',' cluster(s)','  Clean(s)','  visual(s)'])
% disp([round(size(Centers,3)), result, round(sigma), Parameter.values.(Method) ,...
%     Parameter.normalize ,round(itertime), round(cleaningtime),round(visualtime)])

if ischar(Parameter.Context)
    if ischar(Parameter.Spatil.CoOcThr);Parameter.Spatil.CoOcThr=nan;end
    fprintf(strcat('#Centers    psnr    Context    noise    ,',Parameter.Context,'-Lambda  context(s)  K2    NN    CC_{Thr}\n',...
        ' %3u       %2.3f   %2.3f     %3u           %G          %G     %u      %u     %3.2E\n'),...
        size(Centers,3), result,result2, sigma,Parameter.Spatil.lambda,...
        contexttime,length (unique(AssignVec2)),Parameter.Spatil.NN,Parameter.Spatil.CoOcThr);

% disp(['  #Centers    psnr    noise    ',Parameter.Context,'    normalize',' context(s)'])
% disp([length(unique(AssignVec2)), result2, round(sigma), Parameter.Spatil.lambda ,...
%     Parameter.normalize ,round(contexttime)])

end

% rmpath(NewPath);