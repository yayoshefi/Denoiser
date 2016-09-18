function AssignVec2=CoMeans(Data,AssignVec,Centers)
%% AssignVec2=CoMeans(Data,AssignVec,Centers)
% the function return an assignment the same size as AssignVec which is
% based on minimizing both Co-Occurrence entropy and also distance in visul
% space.
% this is consistent with Tokt 1st Scheme of minimizing the functional
global Parameter Analysis

SpatialRefernce='MessagePass';
change=''; cnt=[];   alpha=1.03;

K=size(Centers,3); M=size(Data,2);

iter_step=Analysis.DebuggerIter;
wsize=sqrt(Parameter.wsize2);
lambda=Parameter.spatial.lambda;        rule=Parameter.spatial.UpdateRule;
                                        AssginType=Parameter.spatial.AssginType;
if (Parameter.ORACLE) && (Analysis.ORACLE.level>0)
    Centers=Analysis.ORACLE.Centers;            end
Distances = patch2center (Data,Centers);
[P_i]=affinity (Data, Centers,0,true)';
[CoOc,P_Ni]=lcm(AssignVec,'hard',Parameter.spatial.CoOc);
CoOc_Hold=CoOc_V1(CoOc,false,'entropy');

Lhat=AssignVec;


Parameter.spatial.MaxIter=20;
for iter=1:Parameter.spatial.MaxIter
    if (Parameter.ORACLE) && (Analysis.ORACLE.level>2)
        CoOcThr=Analysis.ORACLE.CoOc;
        [~,P_Ni]=lcm(Lhat,'hard',Parameter.spatial.CoOc);
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

    switch SpatialRefernce
        case 'MessagePass'
            L=(1-lambda)*P_i + lambda*P_Ni*CoOcThr;
        case 'ML'
            L=(1-lambda)*P_i + lambda*P_Ni*CoOcThr';       %notice to invert CC
    end
    Old_Lhat=Lhat;
    [Pr,Lhat]=max(L,[],2);
   
% De-Bugging
    if (iter_step < Parameter.spatial.MaxIter)
    cnt=[cnt,sum(Lhat~=Old_Lhat(:))];
    [CC_Hnew,epsNorm]=CoOc_V1 (CoOcThr,false,'both',Parameter.spatial.CoOcThr);
    ind=sub2ind([K,M],1:M,Lhat);
    score=sum(Distances(ind) )/M+epsNorm;
    Analysis.iterations(iter)=struct('Lhat',Lhat,'CoOc',CoOcThr,...
        'changes',cnt(end),'epsNorm',epsNorm);
    if ~mod(iter,iter_step)              %output to command window
        ratio=CC_Hnew/CoOc_Hold; CoOc_Hold=CC_Hnew;
        fprintf(['iter: %u. %u unique clusters. \tCoOc entropy ratio is %1.3f\n',...
            '%s\t\t\t\t\t\t\t\tCoOc eps norm:%1.4G (i.e. %.3f of the matrix is not 0)\n'],...
            iter,length (unique(Lhat)),ratio,change,epsNorm,(epsNorm/(K^2)));
        fprintf('pixels changed last iterations');disp(cnt)
        change=''; cnt=[];
    
        if Analysis.DebuggerMode 
            Debug(CoOcThr,Lhat,Pr,iter,P_i,P_Ni,L,Analysis.samp);  end
    end
    end

    
%    fixed  Centers-> need to comment next 3 rows speeds up(no calc affinity)
    if rule==1    % update Centers and affinity
        change='centers changed';
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

function[] = Debug (CoOc,Lhat,Pr,iter,S,P_Ni,L,samp)
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

function Distances = patch2center (Data,Centers)
K=size(Centers,3);
Distances=zeros (size(Data,2),K);
for k=1:K
    pointdist=sum(Data.^2)-2*Centers(:,:,k)'*Data+Centers(:,:,k)'*Centers(:,:,k);
    Distances(k,:)=pointdist.^0.5;
end

end