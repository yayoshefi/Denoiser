%%---------------- Script: Pyramid Context-------------------Jan Version
clearvars
load Database;%Images = createImages();
Image=barbara;
%%--------------------------- PARAMETERS ------------------------------
global Parameter Analysis

Method='kmeans';        %Distance  ,  VarianceSplit , kmeans , 'Spectral'
sigma=40;
wsize=9;
normalize=0;            %normalize 0-do nothing ; 1-only bias; 2- bias and gain
metric ='euclidean';    %distance function can be 'euclidean','mahalanobis'
                        %'varing_cluster_size'
Parameter.Spatil.lambda=20;                        

Parameter=struct('row',size(Image,1),'col',size(Image,2),'Method',Method,...
    'sigma',sigma,'wsize2',wsize^2,'normalize',normalize,'metric',metric);

%% ------------------------- INITIALIZATION ---------------------------
NewPath=genpath([pwd,'/Util']);
addpath(NewPath);
setGlobalParameter();

row=Parameter.row;      col=Parameter.col;
Noise=randn(size(Image))*sigma;
Input{1}=double(Image)+Noise;
% ************** Pyramid  *************
maxlevel=2;
for level=2:maxlevel
   Input{level}=impyramid(Input{level-1},'reduce');
end
Parameter.values.kmeans=5;
% ************** Pyramid  *************
d=1;
AssignVec{level}=ones(1,(size(Input{level},1)-wsize+1)*(size(Input{level},1)-wsize+1));
for level=maxlevel:-d:1
    Centers=[];
    if level~=maxlevel; Parameter.values.kmeans=7; end
    Data=im2col(double(Input{level}),[wsize,wsize],'sliding');
    for label=max (AssignVec{level}):-1:1
Data_p=Data(:,AssignVec{level}==label);
X=Data_p;

% X=g*p+b
Xmean=mean(X);
X=X-Xmean(ones(wsize^2,1),:);               % X- bias normalize
Xnorm=(sum(X.^2)).^0.5;
Xn=X./(Xnorm(ones(wsize^2,1),:)+0.001);     % Xn- gain normalize

if Parameter.normalize==2; Patches=Xn;
elseif Parameter.normalize==1; Patches=X; 
else Patches=Data_p;
end
clearvars Xmean Xnorm row col
%% ----------------------- Cluster -----------------------------------
tic
[AssignVec_temp, Centers_temp,Energy,Basis]= FindClusters(Patches,'maxsubspace',Parameter.MSS);
itertime=toc;
AssignVec{level}(AssignVec{level}==label)=(label-1)*Parameter.values.kmeans+AssignVec_temp;
Centers=cat(3,Centers,Centers_temp);
    end
if level~=1
    AssignMat{level}=reshape(AssignVec{level},[sqrt(length(AssignVec{level})),sqrt(length(AssignVec{level}))]);
    AssignMat{level}=imresize(AssignMat{level},size(Input{level-d})-wsize+1,'nearest');
% p=logical(padarray(ones(Parameter.row-wsize+1,Parameter.col-wsize+1),[(wsize-1)/2,(wsize-1)/2]));
    AssignVec{level-d}=AssignMat{level}(:);
else AssignMat{level}=reshape(AssignVec{level},[sqrt(length(AssignVec{level})),sqrt(length(AssignVec{level}))]);
end
end
%% ------------------ Remove Noise --------------------------------
[Output]=removenoise(double(Image),Noise,AssignVec{1});
result=psnr(Output,double(Image),255);
cleaningtime=toc-itertime;

%% -------------------------------------------------------------
%% ------------------  Context  --------------------------------
%% -------------------------------------------------------------
% if ischar(Parameter.Context)
%     [AssignVec2]=SpatialContext (Patches,AssignVec,Centers);
% 
%     [Context_Output]=removenoise(double(Image),Noise,AssignVec2);
% 
%     result2=psnr(Context_Output,double(Image),255);
%     Out={Output,Context_Output};
%     res={result,result2};
% else  AssignVec2=[];
% end
% contexttime=toc-cleaningtime-itertime;

%% ----------------- Visualization ----------------------------
if Analysis.Show
    Analysis.K=size(Centers,3);
    ShowClusters(Image,Noise,AssignVec{1},Centers);
    PrintDnoise ({Output},{result});
    figure;
    for level=1:maxlevel
        subplot(1,maxlevel,level);imagesc(AssignMat{level});
    end

else
PrintDnoise (Out,res,AssignVec{1},AssignVec2)    
end

visualtime=toc-cleaningtime-itertime;


%% ----------------------  Summary  ------------------

disp(['  #Centers    psnr    noise    ',Method,'    normalize',' cluster(s)','  Clean(s)','  visual(s)'])
disp([round(size(Centers,3)), result, round(sigma), Parameter.values.(Method) ,...
    Parameter.normalize ,round(itertime), round(cleaningtime),round(visualtime)])

