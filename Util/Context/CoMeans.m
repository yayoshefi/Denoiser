function AssignVec2=CoMeans(Data,AssignVec,Centers)
%% AssignVec2=CoMeans(Data,AssignVec,Centers)
% the function return an assignment the same size as AssignVec which is
% based on minimizing both Co-Occurrence entropy and also distance in visul
% space.
% this is consistent with Tokt 1st Scheme of minimizing the functional
global Parameter Analysis

SpatialRefernce='MessagePass';      
K=size(Centers,3);                      alpha=1.03;
iter_step=Analysis.DebuggerIter;        wsize=sqrt(Parameter.wsize2);
lambda=Parameter.spatial.lambda;        rule=Parameter.spatial.UpdateRule;
AssginType=Parameter.spatial.AssginType;

%%
if (Parameter.ORACLE) && (Analysis.ORACLE.level>0)
    Centers=Analysis.ORACLE.Centers;            end
[Distances,P_i] = patch2center (Data,Centers,wsize);
[CoOc,P_Ni]=lcm(AssignVec,'hard',Parameter.spatial.CoOc);

Lhat=AssignVec';     Old_Lhat=Lhat;

%% Loop
Parameter.spatial.MaxIter=10;
for iter=1:Parameter.spatial.MaxIter
    if (Parameter.ORACLE) && (Analysis.ORACLE.level>=2)
        CoOcThr=Analysis.ORACLE.CoOc;
        [CoOc,P_Ni]=lcm(Lhat,'hard',Parameter.spatial.CoOc);
    else
        switch AssginType
            case 'soft'
                arg1=reshape(P_i,[],1,K);
            case 'hard'
                arg1=Lhat;
        end
        switch rule
            case {1,2,6}
                [CoOc,P_Ni]=lcm(arg1,AssginType,Parameter.spatial.CoOc);
                CoOcThr=CoOcShrinkage(CoOc);
            case 3
                [CoOc,~]=lcm(arg1,AssginType,Parameter.spatial.CoOc);
                CoOcThr=CoOcShrinkage(CoOc);
            case 4
                 CoOcThr=CoOcShrinkage(CoOc,alpha^iter);
        end
%         if strcmp(AssignType,'soft')
%         [CoOc,P_Ni]=lcm(reshape(P_i,[],1,K),'soft',Parameter.spatial.CoOc);
        
    end
%% De bugging
if Analysis.DebuggerMode 
    Analysis.iterations(iter) = Debugstrct(Lhat,Old_Lhat,CoOc,CoOcThr,...
        Distances,Parameter.spatial.CoOcThr);
    
    if ~mod(iter,iter_step)
        Debugprint (iter)
        if Analysis.Show
            Debug(CoOcThr,Lhat,Pr,iter,P_i,P_Ni,L,Analysis.samp);
        end
    end
end
    
%% Local update    
    switch SpatialRefernce
        case 'MessagePass'
            L=(1-lambda)*P_i + lambda*P_Ni*CoOcThr;
        case 'ML'
            L=(1-lambda)*P_i + lambda*P_Ni*CoOcThr';       %notice to invert CC
    end
    [Pr,Lhat]=max(L,[],2);

%%  Other Rules
    if rule==1    % update Centers and affinity
        [Centers,~,~,~,~]=UpdateCenter(Data,Lhat,false);
       % do I need to update  Lhat in each iter? [Centers,~,Lhat,~,~]
        Centers=cat(3,Centers,inf*ones(wsize^2,1,K-size(Centers,3)));      % To avoid case of degenarated Centers
        [P_i]=affinity (Data, Centers,0,true)';    % maybe: P_i=L
    end
    if rule == 6
        P_i=L;  %TEST 6/22/16 aggregate prob.
    end
end

AssignVec2=Lhat;
end

function CoOcThr=CoOcShrinkage(CoOc,alpha)
global Parameter
% Parameter.spatial.shrink=0.0;
%% zeroize k elements in each row..zs
% element = round (0.4*length(CoOc));
% [SrtCoOc,I]=sort(CoOc);
% ind=sub2ind(size(CoOc),I(:,element)',1:length(CoOc));


per=Parameter.spatial.shrink;
if exist('alpha','var'); per=per*alpha;end

[N,edges] = histcounts(CoOc,10.^(-15:0),'Normalization','cdf');
ind=find (N>per,1);
Parameter.spatial.CoOcThr=edges(ind+1);
%{
if ~isnumeric(Parameter.spatial.CoOcThr)
    switch Parameter.spatial.CoOc
        case 'CC'
            Parameter.spatial.CoOcThr=0.005;
        case 'JP'                           % matrix is normalized row & col
            Parameter.spatial.CoOcThr=1e-9;
        case 'MI'
            Parameter.spatial.CoOcThr=0;     % no shrinkage for Mutual information
        case 'PMI'
             Parameter.spatial.CoOcThr=0;
    end
end
%}
CoOc(CoOc<Parameter.spatial.CoOcThr)=0;
switch Parameter.spatial.CoOc   %Normilazation
    case 'CC'
        CCNorm=sum(CoOc,2);
        CoOcThr=CoOc./CCNorm(:,ones(1,length(CCNorm)),:); CoOcThr(CoOc==0)=0;
    case 'JP'   %matrix is normalized row & col
        CCNorm=sum(CoOc(:));
        CoOcThr=CoOc/CCNorm;
    otherwise   %'MI' or 'PMI'
        CoOcThr=CoOc;   % the mutual information is not a probabilty function
end

end

function [] = Debug (CoOc,Lhat,Pr,iter,S,P_Ni,L,samp)
global Parameter Analysis

ShowProb (cat(3,S,P_Ni*CoOc,L),samp);
subplot(3,1,1);ylabel('Pr. visual');subplot(3,1,2);ylabel('Pr. Context');subplot(3,1,3); xlabel({'pixel',['iter: ',num2str(iter)]})

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
plot(modes);title (strcat('using Thr: ',num2str(Parameter.spatial.CoOcThr)));
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

function [] = Debugprint (iter)
global Analysis
iter_step=Analysis.DebuggerIter;
K=length( Analysis.iterations(iter).CoOc );

cnt  = [Analysis.iterations(iter+1-iter_step:iter).changes];
score= [Analysis.iterations(iter+1-iter_step:iter).score];
l_0=Analysis.iterations(iter).epsNorm;
sparsity= (l_0/(K^2));

fprintf('\niter %i: l_0 norm= %3.0f (%0.3f)\n',iter,l_0,sparsity);
disp( strrep(['pixels changed:' sprintf(' %i,'  ,cnt(:) )], ',', '    '))
disp( strrep(['score  changed:' sprintf(' %3.0f,',score) ], ',', '    '))
end

function [s] = Debugstrct (Lhat,Old_Lhat,CoOc_S,CoOc_T,Distances,Threshold)
K=length(CoOc_S); M=length(Lhat);   beta=0.1;%objective score param

change=sum(Lhat~=Old_Lhat(:));
[CC_Hnew,epsNorm]=CoOc_V1 (CoOc_T,false,'both',Threshold);

ind=sub2ind([K,M],Lhat',1:M);
AvgDist=sum(Distances(ind) )/M;

score=AvgDist+beta*epsNorm;
CoOcConverge = sum( (CoOc_S(:)-CoOc_T(:)).^2 ); %sum( CoOC_T(:) .* ( log(CoOC_T(:)./CoOC_S(:)) ) )
    
s=struct('Lhat',Lhat,'CoOc',CoOc_T,'changes',change,'AvgDist',AvgDist,'epsNorm',epsNorm,...
    'CoOcDist',CoOcConverge,'score',score);
end

function [Distances,E] = patch2center (Data,Centers,wsize)
K=size(Centers,3);
Distances=zeros (K,size(Data,2));
sigma=(wsize^2)/4;
for k=1:K
    pointdist=sum(Data.^2)-2*Centers(:,:,k)'*Data+Centers(:,:,k)'*Centers(:,:,k);
    Distances(k,:)=abs(pointdist).^0.5;
end
E=exp(Distances/(-sigma^2));
cumE=sum(E,1);
NormalizeFactor=1./cumE;
NormalizeFactor( cumE'==0 )=0;
E= bsxfun(@times,E,NormalizeFactor);
E=E';
end