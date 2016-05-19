function [varargout]=ShowCoOc(AssignVec,varargin)
%% function [Output]=ShowCoOc(AssignVec,visual,output)
%
% output  is an optioal field to secifay what will be the return value of the function
% Output = 'CoOc' or 'EpsNorm' 'Entropy'
% 
% visual is an optional bollean argument which specifays wheter to plot the
% Co-Occurence matrix and some more properties.
% in case visual is -[], the default value is true.
%
global Parameter Analysis
row=Parameter.row;      col=Parameter.col;
wsize=sqrt(Parameter.wsize2);
K=Analysis.K;%max(AssignVec)
clstCnt=histcounts(AssignVec,1:K+1);
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

Hprime=DiagonalMult(H,1./HNorm,'l');

Indicator=sparse(AssignVec,1:m*n,ones(1,m*n),K,m*n);
CoOc=Indicator*H;

CoOcNorm=sum(CoOc,2);
CoOcN=CoOc./CoOcNorm(:,ones(1,size(CoOc,1)),:);CoOcN(CoOc==0)=0;

CoOcNprime=diag(1./clstCnt)*Indicator*Hprime;
%% diffrent types of Co-Occurrencce matrix
% h_aH_b=clstCnt'*clstCnt;
% M=CoOc./h_aH_b;                         %like mutual information
% M=CwithWindow (AssignImg)./h_aH_b;
% P=Indicator*Hprime/sum(clstCnt);        %joint probability

%% outputs options
alpha1=0.005; alpha2=0.02; alpha3=0.0001; 
epsNorm1=sum(sum(CoOcN>alpha1));epsNorm2=sum(sum(CoOcN>alpha2)); epsNorm3=sum(sum(CoOcN>alpha3));

LogCoOc=log2(CoOcN);
LogCoOc(CoOcN==0)=0;  % no nan resulting from inf*0;
H_row=-sum( CoOcN.*LogCoOc,2 );
MatrixEntropy=mean(H_row);
%%
visual=true;
if nargin>1
    if ~isempty(varargin{1})
        visual=varargin{1};
    end
end
if nargin>2
    switch varargin{2}
        case 'CoOc'
            varargout{1}=CoOcN;
%             varargout{1}=M;
        case 'EpsNorm'
            varargout{1}=epsNorm1; varargout{2}=epsNorm2;
        case 'Entropy'
            varargout{1}=MatrixEntropy;
    end
end

if visual; Visualize(CoOcN);end  %output figure


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
    xlabel ({strcat( 'amount of labels ',num2str(K) ),strcat( 'active context clusters ',num2str(Analysis.K2) )});axis ([1,K,0,1]); grid minor        
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