function ShowClusters (Image,Noise,AssignVec,Centers,varargin)
%% ShowClusters (Image,Noise,AssignVec,Centers,varargin)
% 
% Image and Noise are same size and when addes they build the data
% if using a Spectral Method with no centers, Use Centers=0


% Parameters
global Parameter Analysis
scrsz = get(groot,'ScreenSize');
figsize=[10 10 1650 1100]; %shold be for 4 images of 512
bins=1000;
NN=Parameter.spatial.NN;  %window in Labels Image

row=Parameter.row;      col=Parameter.col;
wsize=Parameter.wsize2^0.5;
m=row-wsize+1;      n=col-wsize+1;
[VarP,PatchSNR]=ShowPatchSNR (Image,Noise);

p=logical(padarray(ones(m-NN+1,n-NN+1),[floor(NN/2),floor(NN/2)]));
% use p to convert all pathes Images to only inner Image (after NN)

% Initialization
prefix();

if Analysis.Save;  mkdir( strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix) );  end

if ismatrix(Noise)
    Data=im2col(double(Image)+Noise,[wsize,wsize],'sliding');
else
    Data=Image;
end

ClusterMeanVarP=zeros(1,Analysis.K);
ClusterMeanPatchSNR=zeros(1,Analysis.K);

E=false;B=false;
if nargin>4
    Energy=varargin{1};         E=true;
    if nargin>5
        Basis=varargin{2};      B=True;end
end

%% ########################################################################
%% ###################   Lables Analysis          #########################
%% ########################################################################
Analysis.Handles.Lables=figure('Name','Lables analysis ','Position',figsize);%,'outerPosition',scrsz);

if ~verLessThan('matlab','8.4')  % -- Code to run in MATLAB R2014b and later here --
   Analysis.Handles.Lables = (Analysis.Handles.Lables.Number);
end
%% -------- Image cluster labels assignment ----------

subplot(2,2,2);
if length (AssignVec)== (m*n)
    AssignImg=reshape(AssignVec,[m,n]);
else AssignImg=reshape(AssignVec,[m-2,n-2]);
end
imagesc(AssignImg);title ('Image cluster labels assignment');
axis image; %daspect([size(AssignImg) 1])

colormap (Analysis.ColorMap);
%{
%% ---------- Scatter plot in feature space ------------

subplot(2,2,3);

[U,S,~]=svds(Data,3);
Points3D=U(:,1:3)'*Data;

hold on;cm=colormap(jet(K));
scatter(Points3D(1,:),Points3D(2,:),15,AssignVec);
if Centers~=0
    Centers3D=U(:,1:3)'*squeeze(Centers);
    scatter(Centers3D(1,:),Centers3D(2,:),150,'kx');
end
title('scatter plot')
hold off


%% -------- Patches feature space 2D histogrma ----------

subplot(2,2,4)


[H,X,Y]=hist2(Points3D(1,:),Points3D(2,:),500,500);
[H_m,y_i]=max(H);
[~,x_i]=max(H_m);
X1=(X-X(x_i)); Y1=(Y-Y(y_i(x_i)));
imagesc(X,Y,log(H));title('featue space log histogrma plot')

%} 

%% -------- Patches amount for each cluster Histogram ----------------

subplot(2,2,3);
Analysis.density=histc(AssignVec,1:Analysis.K);
bar (Analysis.density); 
hold on;plot(1:Analysis.K,80*ones(1,Analysis.K),'-r');
title('Clusters Histogram')


%% ----------- Energy of each Cluster, Image space --------------
if E
    Energy=sum(Energy,1)/Parameter.wsize2;
    EnergyImg=zeros(size(AssignImg));
end

%% ----------- Patches amount for each cluster, Image spcae --------------
densityImg=zeros(size(AssignImg));

for k=1:Analysis.K
    densityImg(AssignImg==k)=Analysis.density(k);
    if E; EnergyImg(AssignImg==k)=Energy(k); end
    
    ClusterMeanVarP(k)=sum(VarP(AssignVec==k))/Analysis.density(k);
    ClusterMeanPatchSNR(k)=sum(PatchSNR(AssignVec==k))/Analysis.density(k);
     
end
subplot(2,2,1);
imagesc(densityImg);title ('amount of patches for each cluster')
axis image; colorbar

subplot (2,2,4)
range=linspace(min(Analysis.density),max(Analysis.density),50);
ClusterDensity=histc(Analysis.density,range);
bar (range,ClusterDensity,'histc'); 
xlabel ({'Cluster density','# patches belongs to cluster'}); ylabel ('# Clusters')
title('Clusters Histogram')


if Analysis.Save
    saveas(gcf,strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix,...
        '\',num2str(Analysis.Handles.Lables),Analysis.ShortPrefix,'_Lables Analysis.jpg'))
end
if verLessThan('matlab','8.4');     curr=gcf;
else                                curr=get(gcf,'Number'); end
if curr-Analysis.Handles.Lables==Analysis.figures-1;return;end

%% ########################################################################
%% ###################    Clusters properties    ##########################
%% ########################################################################
Analysis.Handles.Clusters=figure ('Name','Clusters properties');%,'outerPosition',scrsz);
if ~verLessThan('matlab','8.4')  % -- Code to run in MATLAB R2014b and later here --
   Analysis.Handles.Clusters = (Analysis.Handles.Clusters.Number);
end
colormap (Analysis.ColorMap);

LabelPatch=im2col(AssignImg,[NN,NN],'sliding');
Analysis.count=zeros(1,size(LabelPatch,2));
Analysis.Neighbours=zeros(1,size(LabelPatch,2));      % in use in figure 3
for k=1:Analysis.K
    Curr=sum(LabelPatch==k);
    Analysis.count=Analysis.count+logical(Curr);
    Analysis.Neighbours(AssignImg(p)==k)=Curr(AssignImg(p)==k);
end
%% ----------- Clusters amount for each mini window --------------
subplot(2,2,1);
CountImg=col2im(Analysis.count,[NN,NN],size(AssignImg),'sliding');

imagesc(CountImg); axis image;title (['amount of diffrent clusters in each ' ,num2str(NN),' window'])
axis off; colorbar

subplot(2,2,2);
NeighboursImg=col2im(Analysis.Neighbours,[NN,NN],size(AssignImg),'sliding');

imagesc(NeighboursImg); axis image;title (['#Similar clusters in each ' ,num2str(NN),' window'])
axis off; colorbar


%% ----------- Energy of each Cluster, Image space --------------

if E 
    subplot(2,2,4)
    imagesc(EnergyImg); axis image;title ('Clusters Energy')
    xlabel({strcat('\Sigma\sigma_{ii}^2 /wsize^2'),'\sigma_{ii}=Var(U(:,i))'})
%     xlabel('sqrt(\Sigma[var_i])');
    axis off; colorbar
end
%%
V_img=zeros (size(AssignImg));     
for k=1:Analysis.K
    ind=AssignVec==k;
    V=sum( sum( (Data(:,ind)-Centers(:,ones(1,sum(ind)),k)).^2,2 )/Parameter.wsize2 )/sum(ind);
    
    V_img(AssignImg==k)=V;
    
end
subplot (2,2,3)
imagesc (V_img), axis image;title ('tmp Energy per cluster, calcualting directly from pixels ')
axis off; colorbar

%{
subplot (3,2,[5,6])
V_vec=sum( Data-Centers(:,AssignVec).^2 )/Parameter.wsize2;
if length (V_vec)== (m*n)
    V_img_pixels=reshape(V_vec,[m,n]);
else V_img_pixels=reshape(V_vec,[m-2,n-2]);
end
imagesc (V_img_pixels), axis image;title ('Energy for each pixel sperattly')
colorbar
%} 
%Per pixel energy need to add another row in subplot

%%%%

%% ----------- Gap distance between NN(1)to NN(2) --------------
if ~Analysis.Fast
    subplot(2,2,3)

    wsize=Parameter.wsize2^0.5;
    [~, Distances]=Dist2SubSpace (Data,Centers,'knn',2); %for simple euclidean
    Gap=Distances(2,:)-Distances(1,:);
    Gap_N=Gap./Distances(1,:);

    Gap_NImg=col2im(Gap_N,[wsize,wsize],[row,col],'sliding');

    imagesc(Gap_NImg);title ('Normalize Gap between the 2 NN of a patch')
    colorbar
    if Parameter.ORACLE
    annotation ('textbox',[0.5,0.1,0.2,0.2],'String','ORACLE','Color','b','FaceAlpha',0.2);end
end
if Analysis.Save
    saveas(gcf,strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix,...
        '\',num2str(Analysis.Handles.Clusters),Analysis.ShortPrefix,'_Cluster Properties_.jpg'))
end
if verLessThan('matlab','8.4');     curr=gcf;
else                                curr=get(gcf,'Number'); end
if curr-Analysis.Handles.Lables==Analysis.figures-1;return;end

%% ########################################################################
%% ###################        Patch Analysis      #########################
%% ########################################################################
Analysis.Handles.Patch=figure ('Name','Patch analysis');%,'outerPosition',scrsz);
if ~verLessThan('matlab','8.4')  % -- Code to run in MATLAB R2014b and later here --
   Analysis.Handles.Patch = (Analysis.Handles.Patch.Number);
end
colormap (Analysis.ColorMap);

%% -------------  Patch Variance (Reference Image)  --------------------
Img_PatchSNR=col2im(PatchSNR,[wsize,wsize],[Parameter.row,Parameter.col]);

Img_VarP=col2im(VarP,[wsize,wsize],[Parameter.row,Parameter.col]);
subplot(2,3,1);
imagesc(Img_VarP); ylabel('patch Variance');  axis off; axis image
colorbar;

%% ---- signal(patch) variance ------------


DataVec=[Img_VarP(p)';Img_PatchSNR(p)';abs(Analysis.DC_Change(p))];
X_DataVec=1:NN^2;
Y= CalcStatistics(DataVec,Analysis.Neighbours,X_DataVec);
Var_range=linspace(min(VarP),max(VarP),bins+1);


[X_Var,Negibours_Var]=ValueVsData(Analysis.Neighbours,Img_VarP(p)',Var_range);


%{
for m=1:NN^2
    E_Var_Signal(m)=mean(Img_VarP(Analysis.Neighbours==m));
    Med_Var_Signal(m)=median(Img_VarP(Analysis.Neighbours==m));
    std_Var_Signal(m)=sqrt(var(Img_VarP(Analysis.Neighbours==m),1));   

    E_PatchSNR(m)=mean(Img_PatchSNR(Analysis.Neighbours==m));
    Med_PatchSNR(m)=median(Img_PatchSNR(Analysis.Neighbours==m));
    std_PatchSNR(m)=sqrt(var(Img_PatchSNR(Analysis.Neighbours==m),1));   
    
    E_DC_Change(m)=mean(DC_Change(Analysis.Neighbours==m));
    Med_DC_Change(m)=median(DC_Change(Analysis.Neighbours==m));
    std_DC_Change(m)=sqrt(var(DC_Change(Analysis.Neighbours==m),1));   
    
end
%}

subplot(2,3,2); title ('analysis for number of patches for each mini window')
plot(X_Var,Negibours_Var)
% plot(X_DataVec,Y.E(1,:),X_DataVec,Y.Med(1,:),X_DataVec,Y.Std(1,:),'cs');
% legend ('Mean','Median','std','Location','NorthEastOutside')

ylabel(['# similar patches in ',num2str(NN),' window']);
xlabel('Var(signal)')

[X_Var,Density_Var]=ValueVsData(densityImg,Img_VarP(p)',Var_range);

% [B,I]=sort(Analysis.density);

subplot(2,3,3); title (' density analysis for each cluster ');
% plot(B,ClusterMeanVarP(I),'r');
plot (X_Var,Density_Var)
ylabel('density _{(#patches in Cluster)} ')
xlabel ('signal variance')


%% ----------- Empirical Patch SNR ---------------

subplot(2,3,4);
imagesc(Img_PatchSNR); ylabel('patch SNR');axis off; axis image
colorbar;

% ---- patch SNR ------------
SNR_range=linspace(min(PatchSNR),max(PatchSNR),bins+1);
[X_SNR,Negibours_SNR]=ValueVsData(Analysis.Neighbours,Img_PatchSNR(p)',SNR_range);
subplot(2,3,5)
plot(X_SNR,Negibours_SNR)
line ([0.45,0.45],[0,max(Negibours_SNR)],'Color',[1,0,0]);

% plot(X_DataVec,Y.E(2,:),X_DataVec,Y.Med(2,:),X_DataVec,Y.Std(2,:),'cs')
% legend ('Mean','Median','std','Location','BestOutside')
ylabel(['# similar patches in ',num2str(NN),' window']);
xlabel('Patch SNR')

subplot(2,3,6);
[X_SNR,Density_SNR]=ValueVsData(densityImg,Img_PatchSNR(p)',SNR_range);

% plot(B,ClusterMeanPatchSNR(I),'r');
plot(X_SNR,Density_SNR)
line ([0.45,0.45],[0,max(Density_SNR)],'Color',[1,0,0]);
ylabel('density _{(#patches in Cluster)} ')
xlabel('[patch SNR');

if Analysis.Save
    saveas(gcf,strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix,...
        '\',num2str(Analysis.Handles.Patch),Analysis.ShortPrefix,'_Plots Analysis_.jpg'))
end
if verLessThan('matlab','8.4');     curr=gcf;
else                                curr=get(gcf,'Number'); end
if curr-Analysis.Handles.Lables==Analysis.figures-1;return;end

%% ########################################################################
%% ###################        Cleaning Analysis      ######################
%% ########################################################################
Analysis.Handles.Cleaning=figure('Name','Cleaning Analysis');%,'outerPosition',scrsz);
if ~verLessThan('matlab','8.4')  % -- Code to run in MATLAB R2014b and later here --
   Analysis.Handles.Cleaning = (Analysis.Handles.Cleaning.Number);
end
colormap (Analysis.ColorMap);

[X_Var,Rmse_Var]=ValueVsData(Analysis.PatchRmse,VarP,linspace(min(VarP),max(VarP),bins+1));

[X_SNR,Rmse_SNR]=ValueVsData(Analysis.PatchRmse,PatchSNR,linspace(min(PatchSNR),max(PatchSNR),bins+1));


%{
UVar=max(VarP); LVar=min(VarP)-0.01;
VarScale=(UVar-LVar)/(bins+1);

USNR=max(PatchSNR); LSNR=min(PatchSNR)-0.01;
SNRscale=(USNR-LSNR)/(bins+1);

X_Var=zeros(1,bins);            X_SNR=zeros(1,bins);
Rmse_Var=zeros(1,bins);         Rmse_SNR=zeros(1,bins);
for i=1:(bins+1)
    X_Var(i)=mean(VarP(  (LVar+(i-1)*VarScale)< VarP & VarP <=(LVar+(i)*VarScale)  ));
    Rmse_Var(i)=mean(PatchRmse(  (LVar+(i-1)*VarScale)< VarP & VarP <=(LVar+(i)*VarScale)  ));
    
    X_SNR(i)=mean(PatchSNR(  (LSNR+(i-1)*SNRscale)< PatchSNR & PatchSNR <=(LSNR+(i)*SNRscale)  ));
    Rmse_SNR(i)=mean(PatchRmse(  (LSNR+(i-1)*SNRscale)< PatchSNR & PatchSNR <=(LSNR+(i)*SNRscale)  ));    
end
%}

C=cov([Analysis.PatchRmse(:),VarP(:),PatchSNR(:)]);
disp (['cov[Rmse,signal Var]=',num2str(C(1,2))]);
disp (['cov[Rmse,Patch SNR]=',num2str(C(1,3))]);

subplot (2,2,1)
plot (X_Var,Rmse_Var)
xlabel('Var(Patch)');           ylabel ('Rmse')

subplot (2,2,3)
plot (X_SNR,Rmse_SNR)
line ([0.45,0.45],[0,max(Rmse_SNR)],'Color',[1,0,0]);
xlabel('Patch SNR');           ylabel ('Rmse')

subplot (2,2,2)
plot(X_DataVec,Y.E(3,:),X_DataVec,Y.Med(3,:),X_DataVec,Y.Std(3,:),'s')
legend ('Mean','Median','std','Location','NorthEastOutside')
ylabel('DC change after Clean]');xlabel(['# similar patches in ',num2str(NN),' window']);

subplot (2,2,4)
imagesc (col2im(abs(Analysis.DC_Change),[wsize,wsize],[row,col],'sliding'))
axis image ;colorbar
title ('DC Change(absolute value) After Clean');

if Analysis.Save
    saveas(gcf,strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix,...
        '\',num2str(Analysis.Handles.Cleaning),Analysis.ShortPrefix,'_Cleaning Analysis_.jpg'))
end
end

% ############### Util Functions   ##############

function ih = PlotLabels(L)

L = single(L);

bL = imdilate( abs( imfilter(L, fspecial('log'), 'symmetric') ) > 0.1, strel('disk', 1));
LL = zeros(size(L),class(L));
LL(bL) = L(bL);
Am = zeros(size(L));
Am(bL) = .5;
ih = imagesc(LL); 
set(ih, 'AlphaData', Am);
colorbar;
colormap 'jet';
end

function Y= CalcStatistics(DataVec,BinBy,X)
n=size (DataVec,1);
E=zeros(n,length(X));       Med=zeros(n,length(X));     Std=zeros(n,length(X));
for x=X
    ind=BinBy==x;
    E(:,x)=mean( DataVec(ind(ones(n,1),:)) );
    Med(x)=median( DataVec(ind(ones(n,1),:) ));
    Std(x)=sqrt(  var( DataVec(ind(ones(n,1),:)),1 )  );
end
Y=struct('E',E,'Med',Med,'Std',Std);
end

function [X,Y]=ValueVsData(Value,Data,binranges)
%% both Value and Data are (1 by pnum)
% X is the mean of each Data point inside the bin

[bincounts,ind]= histc(Data,binranges);
X=zeros(1,length(binranges)-1);         Y=zeros(1,length(binranges)-1);

for i=1:length(binranges)-1
    X(i)= sum( Data(ind==i))/bincounts(i);
    Y(i)= sum(Value(ind==i))/bincounts(i);
    
end
end


function []= prefix()
global Parameter Analysis
% in PrintDnoise there is a subfolder for diffrent noise
Analysis.DirSuffix=strcat('sigma',num2str(Parameter.sigma),'\',Parameter.description,...
    '_Context_',Parameter.Context,...
    '_',num2str(Parameter.spatial.lambda,'%G'),...
    '_CoOc_',Parameter.spatial.CoOc,'_Thr',num2str(Parameter.spatial.CoOcThr,'%G'),...
    Parameter.Method,...
    '_norm',num2str(Parameter.normalize),'_metric_',Parameter.metric);

Analysis.ShortPrefix=strcat('W1-',num2str(sqrt(Parameter.wsize2)),'_Context-',Parameter.Context,...
    '-',num2str(Parameter.spatial.lambda,'%G'),'_NN-',num2str(Parameter.spatial.NN),';');
end
