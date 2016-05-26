function [varargout]=ShowCoOc(AssignVec,visual,output,Method)
%% function [Output]=ShowCoOc(AssignVec,visual,output,Method)
%   Input:
%       visual- is a boolean field that enalbels showing the CC
%       output- 'CoOc'\'EpsNorm'\'Entropy'- chooses what is the output of
%       the function
%       Method- choosed between 
%           'CC' - conditional CC
%           'M'  - Mutual information 
%           'JP' - Joint Probabilty
% 
%   Defaluts:
%       visual is set to true
%       Output is ste to 'CoOc' and method is 'CC' 

global Parameter Analysis
row=Parameter.row;      col=Parameter.col;
wsize=sqrt(Parameter.wsize2);
K=Analysis.K;%max(AssignVec)
if verLessThan('matlab','8.4')
        clstCnt=histc(AssignVec,1:K);
else    clstCnt=histcounts(AssignVec,1:K+1);
end
% clstPr_obj=histogram(AssignVec,1:K+1);
% clstCnt=clstPr_obj.Values;

AssignImg=col2im(AssignVec,[wsize,wsize],[Parameter.row,Parameter.col]);
[m,n]=size(AssignImg);

NN=Parameter.Spatil.NN;  %window2
padding=floor(NN/2);
AssignImg=padarray(AssignImg,[padding,padding],-1);

Neigbour=im2col(AssignImg,[NN,NN],'sliding');
Neigbour(ceil(NN^2/2),:)=[];
H=histc(Neigbour,1:K,1)';
HNorm=sum(H,2);

Indicator=sparse(AssignVec,1:m*n,ones(1,m*n),K,m*n);
CC=Indicator*H;
CCNorm=sum(CC,2);
CCN=CC./CCNorm(:,ones(1,size(CC,1)),:);     CCN(CC==0)=0;

% second way calculating
% Hprime=DiagonalMult(H,1./HNorm,'l');
Hprime=bsxfun(@rdivide,H,HNorm);
Pr_clst=clstCnt/sum(clstCnt);
NullClusters=(clstCnt==0);

CoOcNprime=diag(1./clstCnt)*Indicator*Hprime;   CoOcNprime(NullClusters,:)=0;

%% diffrent types of Co-Occurrencce matrix
h_aH_b=clstCnt'*clstCnt;
M=CC./h_aH_b;                         %like mutual information
% M=CwithWindow (AssignImg)./h_aH_b;
M1=CoOcNprime*diag(1./clstCnt);             M1(:,NullClusters)=0;
P=Indicator*Hprime/sum(clstCnt);        %joint probability
P1=diag(Pr_clst)*CoOcNprime;

if ~exist('Method','var');Method='CC';end
switch Method
    case 'CC'
        CoOc=CoOcNprime;
    case 'M'
        CoOc=M1;
    case 'JP'
        CoOc=P1;    
end

%% outputs options
alpha1=0.005; alpha2=0.02; alpha3=0.0001; 
epsNorm1=sum(sum(CoOc>alpha1));epsNorm2=sum(sum(CoOc>alpha2)); epsNorm3=sum(sum(CoOc>alpha3));

LogCoOc=log2(CoOc);
LogCoOc(CoOc==0)=0;  % no nan resulting from inf*0;
H_row=-sum( CoOc.*LogCoOc,2 );
MatrixEntropy=mean(H_row);
%%
if ~islogical(visual); visual=true; end
if nargout>0
if ~exist('output','var');output='CoOc';end
switch output
    case 'CoOc'
        varargout{1}=CoOc;
    case 'EpsNorm'
        varargout{1}=epsNorm1; varargout{2}=epsNorm2;
    case 'Entropy'
        varargout{1}=MatrixEntropy;
    otherwise
        varargout{1}=CoOc;
end
end

if visual; Visualize(CCN);end  %output figure


    function []=Visualize(CoOcN)
    figure;colormap 'jet'; 
    subplot(2,2,1); imagesc(log(CoOcN+1));
    title ('Co-Occurence matrix');
    xlabel(['mean Entropy per row: ',num2str(MatrixEntropy)],'Color','red')

    subplot(2,2,2);
    modes=sum(CoOcN>0,2);
    plot(modes);title (strcat('using Thr: ','none'));       %num2str(Parameter.Spatil.CoOcThr)));
    grid on; ylabel('||CC_{K,l}||_{\epsilon}')
    xlabel({'Labels K',...
        strcat('\color{blue} \epsilon=',num2str(alpha1),'; |Co-Oc|_{\epsilon} : ',num2str(epsNorm1)),...
        strcat('\color{magenta} \epsilon=',num2str(alpha2),'; |Co-Oc|_{\epsilon} : ',num2str(epsNorm2))});

    subplot(2,2,[3,4]);
    bar( 1:K ,CoOcN(round(K/2),:) )
    title ( strcat('Co-Occurrence prob. for label ', num2str( round(K/2) )) );
    if (isfield(Analysis,'K2')); line2=strcat( 'active context clusters ',num2str(Analysis.K2) );else line2='';end
    xlabel ({strcat( 'amount of labels ',num2str(K) ),line2});axis ([1,K,0,1]); grid minor        
    end
end

function C= CwithWindow (AssignImg)
global Parameter
K=max(AssignImg(:));
NN=Parameter.Spatil.NN;  %window2
padding=(sqrt(Parameter.wsize2)-1)/2;
sigma=1;

C=zeros(K);

Neigbour=im2col(AssignImg,[NN,NN],'sliding')';
[x,y]=meshgrid(-floor(NN/2):floor(NN/2));
Window=exp(-sqrt(x.^2+y.^2)/(sigma^2));
window=fspecial('gaussian',NN);
window=fspecial('average',NN);


window=window(:);   %window(floor(NN^2/2))=[];
L=AssignImg(padding+1:end-padding,padding+1:end-padding);
L=L(:);
for a=1:K
    CenterIsa=Neigbour(L==a,:);
    for b=1:K
        C(a,b)=sum((CenterIsa==b)*window);
    end
end
end

% glcm = graycomatrix(AssignImg,'GrayLimits',[1,K],'NumLevels',K,'Offset',[0,1;1,0],'Symmetric',true);
%CoOc=sum(glcm,3);