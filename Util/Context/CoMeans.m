function AssignVec2=CoMeans(Data,AssignVec,Centers)
%% AssignVec2=CoMeans(Data,AssignVec,Centers)
% the function return an assignment the same size as AssignVec which is
% based on minimizing both Co-Occurrence entropy and also distance in visul
% space.
% this is consistent with Tokt 1st Scheme of minimizing the functional

SpatialRefernce='MessagePass';

global Parameter Analysis
K=size(Centers,3);
wsize=sqrt(Parameter.wsize2);


if (Parameter.ORACLE) && (Analysis.ORACLE.level>0)
        Centers=Analysis.ORACLE.Centers;end
[P_i]=affinity (Data, Centers,0,true)';
Lhat=AssignVec;
NN=Parameter.Spatil.NN;  %window2
padding=floor(NN/2);

ratio=0;            [CC_Hold]=ShowCoOc(Lhat,false,'Entropy',Parameter.Spatil.CoOc);

Parameter.Spatil.MaxIter=30;                    
for iter=1:Parameter.Spatil.MaxIter
    AssignImg=col2im(Lhat,[wsize,wsize],[Parameter.row,Parameter.col]);
    AssignImg=padarray(AssignImg,[padding,padding],-1);

    Neigbour=im2col(AssignImg,[NN,NN],'sliding');
    Neigbour(ceil(NN^2/2),:)=[];
    H=histc(Neigbour,1:K,1)';
    
    if (Parameter.ORACLE) && (Analysis.ORACLE.level>2)
        CCthr=Analysis.ORACLE.CoOc;
    else
        [CC]=ShowCoOc(Lhat,false,'CoOc',Parameter.Spatil.CoOc);
        CCthr=CoOcShrinkage(CC);
        
        [CC_Hnew,epsNorm]=CC_Entropy(CCthr); ratio=CC_Hnew/CC_Hold; CC_Hold=CC_Hnew;
        if ~mod(iter-1,7)
        fprintf(['iter: %u. %u unique clusters. \tCoOc entropy ratio is %1.3f\n',...
            '\t\t\t\t\t\t\t\tCoOc eps norm:%1.4G (i.e. %.3f of the matrix is not 0)\n'],...
            iter,length (unique(Lhat)),ratio,epsNorm,(epsNorm/(K^2)));
        end
    end
%     [I, PixelTotalI]=tf_idf (H);
%     PixelTotalI=(PixelTotalI<100);
%     PixelWeightI=PixelTotalI(:,ones(1,K));
    switch SpatialRefernce
        case 'MessagePass'
            HNorm=sum(H,2);
            E_h=DiagonalMult(H,1./HNorm,'l'); %for partial histograms
%             E_h=E_h.*PixelWeightI;%%
            L=(1-Parameter.Spatil.lambda)*P_i + Parameter.Spatil.lambda*E_h*CCthr;
        case 'CenterPixel'
            E_h=HistDist (H,CCthr,Centers);
            L=(1-Parameter.Spatil.lambda)*P_i + Parameter.Spatil.lambda*E_h*CCthr;
        case 'ML'
            HNorm=sum(H,2);
            E_h=DiagonalMult(H,1./HNorm,'l'); %for partial histograms
            L=(1-Parameter.Spatil.lambda)*P_i + Parameter.Spatil.lambda*E_h*CCthr';       %notice to invert CC
    end
    [Pr,Lhat]=max(L,[],2);
   
    % De-Bugging
    if Analysis.DebuggerMode && ~(mod(iter-1,10));
        Debug(CCthr,Lhat,Pr,iter,P_i,E_h,L,Analysis.samp);end
    Analysis.iterations(iter).AssignVec2=single(Lhat);
    Analysis.iterations(iter).CoOc=single(CCthr);

    % update Centers and affinity
%    fixed  Centers-> need to comment next 3 rows speeds up(no calc affinity)
    if ~mod(iter-1,5)
        [Centers,~,~,~,~]=UpdateCenter(Data,Lhat,false);
       % do I need to update  Lhat in each iter? [Centers,~,Lhat,~,~]
        Centers=cat(3,Centers,inf*ones(wsize^2,1,K-size(Centers,3)));      % To avoid case of degenarated Centers
        [P_i]=affinity (Data, Centers,0,true)';    % maybe: P_i=L
    end
    P_i=L;  %TEST 6/22/16 aggregate prob.
end

AssignVec2=Lhat;
end

function CoOcThr=CoOcShrinkage(CoOc)
global Parameter
if ~isnumeric(Parameter.Spatil.CoOcThr)
    switch Parameter.Spatil.CoOc
        case 'CC'
            Parameter.Spatil.CoOcThr=0.005;
        case 'JP'                           % matrix is normalized row & col
            Parameter.Spatil.CoOcThr=1e-9;
        case 'M'
            Parameter.Spatil.CoOcThr=0;     % no shrinkage for Mutual information
    end
end
CoOc(CoOc<Parameter.Spatil.CoOcThr)=0;
switch Parameter.Spatil.CoOc
    case 'CC'
        CCNorm=sum(CoOc,2);
        CoOcThr=CoOc./CCNorm(:,ones(1,length(CCNorm)),:); CoOcThr(CoOc==0)=0;
    case 'JP'   %matrix is normalized row & col
        CCNorm=sum(CoOc(:));
        CoOcThr=CoOc/CCNorm;
    otherwise   %'M'
        CoOcThr=CoOc;   % the mutual information is no a probabilty function
end

end
function E_h=HistDist (H,CC,Centers)
global Parameter

pnum=size(H,1);         speed='fast';

if Parameter.wsize2==1;     Centers=squeeze(Centers)';
else                        Centers=squeeze(Centers);   end
PrH=DiagonalMult(H,1./sum(H,2),'l');
PrH(H==0)=0;

D_h=zeros(size(H)); E_h=zeros(size(H));
switch speed
    case 'normal'
        for k=1:length (CC)
            for i=1:pnum
                [~,D_h(i,k)]=emd(Centers',Centers',PrH(i,:)',CC(k,:)');
            end
        end
    case 'fast'
        dist=FastEMD(Centers,CC(1,:),Centers,PrH);
        E_h(:,1)=1;
        for k=2:length (CC)
%             D_h(:,k)=FastEMD(Centers,CC(k,:),Centers,PrH); with
%             k=1:length
            tempdist=FastEMD(Centers,CC(k,:),Centers,PrH);
            E_h(tempdist<dist,k)=1;
            E_h(tempdist<dist,1:k-1)=0;
            dist(tempdist<dist)=tempdist(tempdist<dist);
        end
end
        

% affinity=exp(-D_h/(2*Parameter.Spatil.sigma^2));
% E_h=DiagonalMult(affinity,1./sum(affinity,2),'l');
% E_h(affinity==0)=0;

end

function[] = Debug (CoOc,Lhat,Pr,iter,S,E_h,L,samp)
global Parameter Analysis

ShowProb (cat(3,S,E_h*CoOc,L),samp);
subplot(3,1,1);ylabel('Pr. visual');subplot(3,1,2);ylabel('Pr. Hist');subplot(3,1,3); xlabel({'pixel',['iter: ',num2str(iter)]})

wsize=sqrt(Parameter.wsize2);

% ShowCoOc(Lhat); set(gcf,'Name',strcat('Co-Occurence ',num2str(iter), ' iteration'));
Pr_Img=col2im(Pr,[wsize,wsize],[Parameter.row,Parameter.col]);
tmp_Labels=col2im(Lhat,[wsize,wsize],[Parameter.row,Parameter.col]);
if isfield (Analysis,'ORACLE')
    RI=RandIndex(Analysis.ORACLE.AssignVec(:),Lhat(:));
end
figure('Name',strcat('temporal image properties ',num2str(iter), 'iteration'));
subplot(2,2,3);
imshow(CoOc,[]);colormap (Analysis.ColorMap); title ('Co-Occurence matrix');

LCoOc=log2(CoOc);
LCoOc(CoOc==0)=0;  % no nan resulting from inf*0;
H_row=-sum( CoOc.*LCoOc,2 );
H=mean(H_row);

xlabel(strcat('mean entropy for each row is:  ',num2str (H)));
subplot(2,2,4);
modes=sum(CoOc>0,2);
plot(modes);title (strcat('using Thr: ',num2str(Parameter.Spatil.CoOcThr)));
axis([1,size(CoOc,1),0,15]);grid on
xlabel('Labels');ylabel('||_{\epsilon} per Conditonal Pr.')

subplot(2,2,1);
imshow(Pr_Img,[]);title ('Probabilty to be in Lhat'); colormap jet
DrawPixels (gcf, samp,false)
xlabel(strcat('mean prob.= ',num2str(mean(Pr)) ));colorbar
subplot(2,2,2)
imshow(tmp_Labels,[]);title (['Labels iter: ',num2str(iter)]); colormap jet
DrawPixels (gcf, samp,false)
if isfield (Analysis,'ORACLE')
    xlabel({strcat(num2str(length(unique(Lhat))) ,' diffrent labels'),...
        strcat('Rand index is: ',num2str(RI))})
else
    xlabel(strcat(num2str(length(unique(Lhat))) ,' diffrent labels'));
end
end

function [H,varargout]=CC_Entropy(CoOcN,K,m,n)
global Parameter

if isvector(CoOcN)         %case this is AssignVec and not CoOc
    wsize=sqrt(Parameter.wsize2); NN=Parameter.Spatil.NN; padding=floor(NN/2);
    AssignImg=col2im(CoOcN,[wsize,wsize],[Parameter.row,Parameter.col]);
    AssignImg=padarray(AssignImg,[padding,padding],-1);

    Neigbour=im2col(AssignImg,[NN,NN],'sliding');
    Neigbour(ceil(NN^2/2),:)=[];
    H=histc(Neigbour,1:K,1)';
    
    Indicator=sparse(CoOcN,1:m*n,ones(1,m*n),K,m*n);
    CC=Indicator*H;
    CCNorm=sum(CC,2);
    CoOcN=CC./CCNorm(:,ones(1,K),:); CoOcN(CC==0)=0; %to avoid 0/0=nan
end

LogCoOc=log2(CoOcN);
LogCoOc(CoOcN==0)=0;  % no nan resulting from inf*0;
H_row=-sum( CoOcN.*LogCoOc,2 );
H=mean(H_row);
if nargout>1
     varargout{1}=sum(sum(CoOcN>Parameter.Spatil.CoOcThr));
end
end