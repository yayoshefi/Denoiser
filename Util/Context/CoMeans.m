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

if Parameter.ORACLE
        Centers=Analysis.ORACLE.Centers;end
[S]=affinity (Data, Centers,0,true)';
Lhat=AssignVec;
NN=Parameter.Spatil.NN;  %window2
padding=floor(NN/2);
%Analysis.LabelsSize=Analysis.LabelsSize-(NN-1);
m=Analysis.LabelsSize(1)    ;n=Analysis.LabelsSize(2);
%p=logical(padarray(ones(m,n),[padding,padding]));
%         
%Analysis.LabelsSize=Analysis.LabelsSize+(NN-1); %restore values to origin

samp=randperm(size(S,1),3);
for iter=1:4
    
%     if iter>1; SpatialRefernce='emd';end
    AssignImg=col2im(Lhat,[wsize,wsize],[Parameter.row,Parameter.col]);
    AssignImg=padarray(AssignImg,[padding,padding],-1);

    Neigbour=im2col(AssignImg,[NN,NN],'sliding');
    Neigbour(ceil(NN^2/2),:)=[];
    H=histc(Neigbour,1:K,1)';
    
    if Parameter.ORACLE
        CCthr=Analysis.ORACLE.CoOc;
    else
        if ~isinteger(Parameter.Spatil.CoOcThr);Parameter.Spatil.CoOcThr=0.005;end
        Indicator=sparse(Lhat,1:m*n,ones(1,m*n),K,m*n);
        CC=Indicator*H;
        CCNorm=sum(CC,2);
        CCN=CC./CCNorm(:,ones(1,K),:); CCN(CC==0)=0; %to avoid 0/0=nan

        CCN(CCN<Parameter.Spatil.CoOcThr)=0;
        CCNorm=sum(CCN,2);
        CCthr=CCN./CCNorm(:,ones(1,K),:); CCthr(CCN==0)=0;
    end
    
    switch SpatialRefernce
        case 'MessagePass'
            E_h=H/(NN^2-1);
        case 'emd'
            E_h=HistDist (H,CC,Centers);
    end
            
    L=(1-Parameter.Spatil.lambda)*S + Parameter.Spatil.lambda*E_h*CCthr;
    % when not caclculating partial histogrmas
    %L=S;
    %L(p,:)=(1-Parameter.Spatil.lambda)*reshape(S(p,:),m*n,K)+Parameter.Spatil.lambda/(NN^2-1)*H*CCthr;
    
    [Pr,Lhat]=max(L,[],2);
    if Analysis.DebuggerMode && ~(mod(iter+1,2))
        ShowProb (cat(3,S,E_h,L),samp);Debug(CCthr,Lhat,Pr,iter);end

%    fixed  Centers
   [Centers,~,Lhat,~,~]=UpdateCenter(Data,Lhat,false);
   [S]=affinity (Data, Centers,0,true)';
end

AssignVec2=Lhat;
end
function E_h=HistDist (H,CC,Centers)
global Parameter

pnum=size(H,1);         speed='fast';

if Parameter.wsize2==1;     Centers=squeeze(Centers)';
else                        Centers=squeeze(Centers);   end
PrH=DiagonalMult(H,1./sum(H,2),'l');
PrH(H==0)=0;

D_h=zeros(size(H));
switch speed
    case 'normal'
        for k=1:length (CC)
            for i=1:pnum
                [~,D_h(i,k)]=emd(Centers',Centers',PrH(i,:)',CC(k,:)');
            end
        end
    case 'fast'
        for k=1:length (CC)
            D_h(:,k)=FastEMD(Centers,CC(k,:),Centers,PrH);
        end
end
        

affinity=exp(-D_h/(2*Parameter.Spatil.sigma^2));
E_h=DiagonalMult(affinity,1./sum(affinity,2),'l');
E_h(affinity==0)=0;

end

function[] = Debug (CoOc,Lhat,Pr,iter)
global Parameter
wsize=sqrt(Parameter.wsize2);

% ShowCoOc(Lhat); set(gcf,'Name',strcat('Co-Occurence ',num2str(iter), ' iteration'));
Pr_Img=col2im(Pr,[wsize,wsize],[Parameter.row,Parameter.col]);
tmp_Labels=col2im(Lhat,[wsize,wsize],[Parameter.row,Parameter.col]);

figure('Name',strcat('temporal image properties ',num2str(iter), 'iteration'));
subplot(2,2,3);
imagesc(log(CoOc+1));colormap jet; title ('Co-Occurence matrix');

LCoOc=log2(CoOc);
LCoOc(CoOc==0)=0;  % no nan resulting from inf*0;
H_row=-sum( CoOc.*LCoOc,2 );
H=mean(H_row);

xlabel(strcat('mean entropy for each row is:  ',num2str (H)));
subplot(2,2,4);
modes=sum(CoOc>0,2);
plot(modes);title (strcat('using Thr: ',num2str(Parameter.Spatil.CoOcThr)));
axis([1,size(CoOc,1),0,15]);grid on
xlabel('Labels');ylabel('modes for every cluster')

subplot(2,2,1);
imagesc(Pr_Img);title ('Probabilty to be in Lhat'); colormap jet
xlabel(strcat('mean prob.= ',num2str(mean(Pr)) ));colorbar
subplot(2,2,2)
imagesc(tmp_Labels);title ('Temp Label image'); colormap jet
xlabel(strcat(num2str(length(unique(Lhat))) ,' diffrent labels'))

end
