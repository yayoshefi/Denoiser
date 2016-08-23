% %%---------------- Data Collecting-------------------Git Version
clearvars ContextPsnr ClusterCompre Psnr  Labeling RulesPsnr RulesClusterCompre full_Data
load ExpImages.mat;     load rectImage.mat

description='kmeans compare2 gabor';
%% --------------------------- PARAMETERS ------------------------------
global Parameter Analysis
%check struct Name
Method='kmeans'; 
metric ='euclidean';

sigma=50;
wsize=11;       NN=9;

% ---- arrays ----
lam=3*10.^(-5:-1);
% lambda_array=[0.001,];
lambda_array=[lam,flip(1-lam)];
rules_array=[3,4];
CoOcType_array={'CC','MI'};
AssignType_array={'hard'};
ShrinkPer=[0,0.025,0.15];%[0.0005,0.0001];
L=length(lambda_array);     T=length(CoOcType_array);

Arrays=struct('Lambda',lambda_array,'CoOcPerThr',ShrinkPer);

for i=14:14
    Image=Irect(i).Image;
    name=(Irect(i).name);   disp(name)
Parameter=struct('description',description,'ImageName',name,'row',size(Image,1),'col',...
    size(Image,2),'Method',Method,'sigma',sigma,'wsize2',wsize^2,'normalize',...
    0,'metric',metric);
%% ------------------------- INITIALIZATION ---------------------------
NewPath=genpath([pwd,'/Util']);
addpath(NewPath);
setGlobalParameter();                   Analysis.Show=true;
row=Parameter.row;                      col=Parameter.col;
Parameter.spatial.NN=NN;
Analysis.DebuggerIter=1000;
ImageDir=strcat(Parameter.location,'\Results\',date,'\',Parameter.ImageName,'_sigma',num2str(Parameter.sigma),'\');


Analysis.Arrays=Arrays;
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
% DM3D
[BM3dresult, BM3dOutput] = BM3D(im2double(Image), im2double(Image)+(Noise/255), sigma,'np',0);

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
tic
[AssignVec, Centers,Energy,Basis]=...
    FindClusters(Patches,'maxsubspace',Parameter.MSS);
itertime=toc;
            
%% ------------------ Remove Noise --------------------------------

[Output]=removenoise(double(Image),Noise,AssignVec);
result=psnr(Output,double(Image),255);
cleaningtime=toc-itertime;
            
%% ------------------  Context  --------------------------------
best='';
for t=1:length(CoOcType_array);
    Parameter.spatial.CoOc=CoOcType_array{t};
    best=[best,CoOcType_array{t},'\n'];
for rule=rules_array;
    Parameter.spatial.UpdateRule=rule;      Parameter.description=['rule #',num2str(rule)];
                                            best=[best,'rule #',num2str(rule),':'];
for a=1:length(AssignType_array);
    Parameter.spatial.AssginType=AssignType_array{a};
    disp(['rule ',num2str(rule),', ',CoOcType_array{t},', ',AssignType_array{a},' assignment'])
    clearvars ContextPsnr ClusterCompre
    
    
for thr=1:length(ShrinkPer);
    Parameter.spatial.shrink=ShrinkPer(thr);
for l=1:length(lambda_array)
    Parameter.spatial.lambda=lambda_array(l);
    
    startcontext=toc;
    [AssignVec2]=SpatialContext (Patches,AssignVec,Centers);
    contexttime=toc-startcontext;
    
    samp_iter=find (( [Analysis.iterations.changes]/length(AssignVec2) )<0.002,1) ;
    if isempty(samp_iter); samp_iter=5;end
    % find the iteration where less the 0.2 % of the image changed labels
    
    [Context_Output]     =removenoise(double(Image),Noise,AssignVec2);
    [samp_iter_Output   ]=removenoise(double(Image),Noise,Analysis.iterations(samp_iter).Lhat);
    
    result2      =psnr(Context_Output,double(Image),255);
    samp_iter_result2=psnr(samp_iter_Output   ,double(Image),255);
    
    [RI,MH  ]           = RandIndex(AssignVec (:),ORACLE(:));
    [RI2,MH2]           = RandIndex(AssignVec2 (:),ORACLE(:));
    [samp_RI,samp_MH]   = RandIndex(Analysis.iterations(samp_iter).Lhat (:),ORACLE(:));
%% structures to save
    temppsnr=struct('lambda',lambda_array(l),'CoOcPerThr',Parameter.spatial.shrink,...
        'MaxIter',Parameter.spatial.MaxIter,'psnr2',result2,'dev2',result2-result,...
        'samp_iter',samp_iter,'samp_result',samp_iter_result2,'samp_dev',samp_iter_result2-result,'Thr',Parameter.spatial.CoOcThr);
    tempCmpClst=struct('lambda',lambda_array(l),'CoOcPerThr',Parameter.spatial.shrink,...
        'MH',MH,'MH2',MH2,'samp_iter',samp_iter,'samp_MH',samp_MH);
    
    ContextPsnr(l+(thr-1)*L)=temppsnr;
    
    ClusterCompre(l+(thr-1)*L)=tempCmpClst;
    
    Labeling(l+(thr-1)*L)=struct('lambda',lambda_array(l),'CoOcPerThr',ShrinkPer(thr),'AssignVec2',AssignVec2);

%% Summary
    fprintf(['#Centers    psnr    Context    noise    ',Parameter.Context,'-Lambda  context(m)  K2   NN   CoOc    CC_{Per}   psnr@%i\n',...
        ' %3u       %2.3f   %2.3f     %3u            %3G         %2.2G   %u      %u    %s     %2.2G         %2.3f\n'],...
        samp_iter,round(size(Centers,3)), result,result2, sigma,Parameter.spatial.lambda,...
        contexttime/60,length (unique(AssignVec2)),NN,Parameter.spatial.CoOc,Parameter.spatial.shrink,samp_iter_result2)

    fprintf (['Basic        @%i            @%i             ORACLE       BM3D\n',...
                '%2.3f      %1.4f       %1.4f       %2.3f (%1.4f)   %2.3f \n\n'],...
                samp_iter,Parameter.spatial.MaxIter,result,...
                samp_iter_result2-result,result2-result,ORACLEresult,ORACLEresult-result,BM3dresult)
    close all
    PrintDnoise ({Output,Context_Output},{result,result2,ORACLEresult},AssignVec,AssignVec2)
    
end % Lambda
end % CoOcThr
[value2,index2]         =max([ContextPsnr.psnr2]);
[samp_value,samp_index] =max([ContextPsnr.samp_result]);
if (samp_index~=index2); disp(['mismatch at rule: ',num2str(rule)]);end
if samp_value<value2
        ind=index2;
else    ind=samp_index;
end
CoOc_V1 (lcm(Labeling(ind).AssignVec2,Parameter.spatial.AssginType,Parameter.spatial.CoOc),true)
saveas (gcf,strcat(ImageDir,'CoOc.png'))

ContextPsnr(ind).CoOc=CoOcType_array{t};        ClusterCompre(ind).CoOc=CoOcType_array{t};
RulesPsnr(rule+(t-1)*T)=ContextPsnr(ind);
RulesClusterCompre(rule+(t-1)*T)=ClusterCompre(ind);

Save4Latex('name',[name,'_'],'lambda',ContextPsnr(ind).lambda,'CoOc',ContextPsnr(ind).CoOc,...
    'Thr',ContextPsnr(ind).Thr,'rule',rule,'sigma',sigma,'Method',Method)
best=strcat(best,'lam=',num2str(ClusterCompre(ind).lambda,'%1.3E') );
end % AssignType
end % rules
SaveStaticImg(ORACLE,ORACLEOutput,ORACLEresult,BM3dresult, BM3dOutput,...
                name,sigma,{Image,Input},{'Original', 'Noisy'})
end % CoOcType

Ref(1)=struct('wsize',wsize,'sigma',sigma,'Psnr',result,'Method',Method,'AssignVec',[]);
Ref(2)=struct('wsize',wsize,'sigma',sigma,'Psnr',ORACLEresult,'Method','ORACLE','AssignVec',ORACLE);
        % ORACLE
        Labeling(thr*L+1).AssignVec2=ORACLE;
        Labeling(thr*L+1).lambda='ORACLE';
              
        full_Data(i)=struct('ImageName',name,'Reference',Ref,...
            'ContextPsnr',RulesPsnr,'ClusterCompre',RulesClusterCompre,'Arrays',Arrays);
end % Image
 %% save info
save (strcat(Parameter.location,'\Results/',date,'/','full_Data.mat'), 'full_Data', '-v7.3')
frprintf(best)
% for i=1:size(full_Data,2)
%     disp([full_Data(i).ImageName,'  ','rule #',num2str(rule)])
%     disp(full_Data(i).ContextPsnr(rule))
% end