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
AssginType=Parameter.CoOc.AssginType;

%%
if (Parameter.ORACLE) && (Analysis.ORACLE.level>0)
    Centers=Analysis.ORACLE.Centers;            end
[Distances,P_i] = patch2center (Data,Centers,wsize);
[CoOc_S,P_Ni]=lcm(AssignVec,'hard',Parameter.CoOc.Type);

Lhat=AssignVec';     Old_Lhat=Lhat;

d_0= sum(Distances( sub2ind(size(Distances),Lhat',1:length(Lhat)) ) )/length(Lhat);
d_iter=d_0;

%% Loop
Parameter.spatial.MaxIter=10;               iter=0;
% while (iter < Parameter.spatial.MaxIter) && ( abs(d_iter - d_0) < 0.1 ) 
% iter=iter+1;
for iter=1:Parameter.spatial.MaxIter
%% evolving lambda
    if ( ~mod(iter-1,3) );lambda = lambda*10;
    [~,P_Ni]=lcm(Lhat,'hard',Parameter.CoOc.Type);
% %     Parameter.CoOc.ShrinkPer=Parameter.CoOc.ShrinkPer+0.075;
    end
%% evolving lambda

    if (Parameter.ORACLE) && (Analysis.ORACLE.level>=2)
        CoOc_T=Analysis.ORACLE.CoOc;
        [CoOc_S,P_Ni]=lcm(Lhat,'hard',Parameter.CoOc.Type);
    else
        switch AssginType
            case 'soft'
                arg1=reshape(P_i,[],1,K);
            case 'hard'
                arg1=Lhat;
        end
        switch rule
            case {1,2,6}
                [CoOc_S,P_Ni]=lcm(arg1,AssginType,Parameter.CoOc.Type);
                CoOc_T=CoOcShrinkage(CoOc_S);
            case 3
                [CoOc_S,~]=lcm(arg1,AssginType,Parameter.CoOc.Type);
                CoOc_T=CoOcShrinkage(CoOc_S);
            case 4
                 CoOc_T=CoOcShrinkage(CoOc_S,alpha^iter);
        end
%         if strcmp(AssignType,'soft')
%         [CoOc,P_Ni]=lcm(reshape(P_i,[],1,K),'soft',Parameter.CoOc.Type);
        
    end
%% De bugging
d_iter= sum(Distances( sub2ind(size(Distances),Lhat',1:length(Lhat)) ) )/length(Lhat);

if Analysis.DebuggerMode 
    Analysis.iterations(iter) = Debugstrct(Lhat,Old_Lhat,CoOc_S,CoOc_T,...
        Distances,Parameter.CoOc.epsilon,Analysis.LabelsSize);
    
    if ~mod(iter,iter_step)
        Debugprint (iter)
        if Analysis.Show
            DebugShow(CoOc_T,Lhat,Pr,iter,P_i,P_Ni,L,Analysis.samp);
        end
    end
end
    
%% Local update    
    switch SpatialRefernce
        case 'MessagePass'
            L=(1-lambda)*P_i + lambda*P_Ni*CoOc_T;
        case 'ML'
            L=(1-lambda)*P_i + lambda*P_Ni*CoOc_T';       %notice to invert CC
    end
    Old_Lhat=Lhat;
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
if Analysis.DebuggerMode ;
    Analysis.iterations(end+1) = Debugstrct(Lhat,Old_Lhat,lcm(Lhat,AssginType,Parameter.CoOc.Type),CoOc_T,...
        Distances,Parameter.CoOc.epsilon,Analysis.LabelsSize);
end
AssignVec2=Lhat;
end

function CoOcThr=CoOcShrinkage(CoOc,alpha)
global Parameter
% Parameter.CoOc.ShrinkPer=0.0;
ShrinkageType  = Parameter.CoOc.ShrinkType;

switch ShrinkageType
    case 'row'
        %% Zeroize k elements in each row;
        k = round (Parameter.CoOc.ShrinkPer * length(CoOc));
        if k>0
            [SrtCoOc]=sort(CoOc,2);
            CoOc ( bsxfun(@ge,CoOc,SrtCoOc(:,k)) )=0;
        end
        
    case 'matrix'
        %% Zeroize bottom precentage of elemnts (entire matrix)
        per=Parameter.CoOc.ShrinkPer;
        if exist('alpha','var'); per=per*alpha;end
        switch Parameter.CoOc.Type  
            case 'CC';                bins=10.^(-5:(1/3):0);
            case 'MI';                bins=10.^(-15:0);
            otherwise;                bins=10.^(-15:0);
        end
        [N,edges] = histcounts(CoOc,bins,'Normalization','cdf');
        ind=find (N>per,1);
        Parameter.CoOc.Thr=edges(ind+1);
 
        CoOc(CoOc<Parameter.CoOc.Thr)=0;
    case 'none'
    case 'epsilon'
        %% Zeroize by constant factor epsilon
        Parameter.CoOc.Thr=Parameter.CoOc.epsilon;
        CoOc(CoOc<Parameter.CoOc.Thr)=0;
end

%%  Normilazation
switch Parameter.CoOc.Type  
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

function [] = DebugShow (CoOc,Lhat,Pr,iter,S,P_Ni,L,samp)
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
plot(modes);title (strcat('using Thr: ',num2str(Parameter.CoOc.Thr)));
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
l_0=Analysis.iterations(iter).l_0;
sparsity= (l_0/(K^2));

fprintf('\niter %i: l_0 norm= %3.0f (%0.3f)\n',iter,l_0,sparsity);
disp( strrep(['pixels changed:' sprintf(' %i,'  ,cnt(:) )], ',', '    '))
disp( strrep(['score  changed:' sprintf(' %3.0f,',score) ], ',', '    '))
end

function [s] = Debugstrct (Lhat,Old_Lhat,CoOc_S,CoOc_T,Distances,epsilon,rwcl)
K=length(CoOc_S); M=length(Lhat);   beta=0.1;%objective score param

img=reshape( uint16(Lhat) ,rwcl);
change=sum(Lhat~=Old_Lhat(:));
[CC_Hnew,l_0,l_1]=CoOc_V1 (CoOc_S,false,'both',epsilon);

ind=sub2ind([K,M],Lhat',1:M);
AvgDist=sum(Distances(ind) )/M;

score=AvgDist+beta*l_1;
CoOcConverge = sum( (CoOc_S(:)-CoOc_T(:)).^2 );
d1 = ( log(CoOc_T(:)./CoOc_S(:)) );         d1( isnan(d1) )=0;      d1( isinf(d1) )=0; %contradiction to divergance definition
d2 = ( CoOc_T(:) .* d1 );                   d2 ( isnan(d2) )=0;
divergance = sum( d2 );
    
s=struct('Lhat',img,'CoOc',CoOc_T,'changes',change,'AvgDist',AvgDist,'l_0',l_0,...
    'l_1',l_1,'CoOcDist',CoOcConverge,'divergance',divergance,'score',score);
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