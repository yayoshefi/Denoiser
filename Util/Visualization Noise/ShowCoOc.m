function [varargout]=ShowCoOc(AssignVec)
global Parameter
row=Parameter.row;      col=Parameter.col;
wsize=sqrt(Parameter.wsize2);
K=max(AssignVec);
alpha1=0.005; alpha2=0.02;

AssignImg=col2im(AssignVec,[wsize,wsize],[Parameter.row,Parameter.col]);
[m,n]=size(AssignImg);

% NN=Parameter.Spatil.NN;
NN=5;
padding=floor(NN/2);
AssignImg=padarray(AssignImg,[padding,padding],-1);

Neigbour=im2col(AssignImg,[NN,NN],'sliding');
Neigbour(ceil(NN^2/2),:)=[];
H=histc(Neigbour,1:K,1)';

Indicator=sparse(AssignVec,1:m*n,ones(1,m*n),K,m*n);
CoOc=Indicator*H;


CoOcNorm=sum(CoOc,2);
CoOcN=CoOc./CoOcNorm(:,ones(1,size(CoOc,1)),:);CoOcN(CoOc==0)=0;

LogCoOc=log2(CoOcN);
LogCoOc(CoOcN==0)=0;  % no nan resulting from inf*0;
H_row=-sum( CoOcN.*LogCoOc,2 );
H=mean(H_row);

epsNorm1=sum(sum(CoOcN>alpha1));epsNorm2=sum(sum(CoOcN>alpha2));

figure;colormap 'jet'; 
subplot(1,2,1); imagesc(log(CoOcN+1));
ylabel ('Co-Occurence matrix');
xlabel(['mean Entropy per row: ',num2str(H)],'Color','red')

subplot(1,2,2);
modes=sum(CoOcN>0,2);
plot(modes);title (strcat('using Thr: ',num2str(Parameter.Spatil.CoOcThr)));
grid on; ylabel('modes for every cluster')
xlabel({'Labels',...
    strcat('\color{blue} \epsilon=',num2str(alpha1),'; |Co-Oc|_{\epsilon} : ',num2str(epsNorm1)),...
    strcat('\color{blue} \epsilon=',num2str(alpha2),'; |Co-Oc|_{\epsilon} : ',num2str(epsNorm2))});

varargout{1}=epsNorm1; varargout{2}=epsNorm2;
end

% glcm = graycomatrix(AssignImg,'GrayLimits',[1,K],'NumLevels',K,'Offset',[0,1;1,0],'Symmetric',true);
%CoOc=sum(glcm,3);