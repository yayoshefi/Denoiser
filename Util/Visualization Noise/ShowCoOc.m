function ShowCoOc(AssignVec)
global Parameter
row=Parameter.row;      col=Parameter.col;
wsize=sqrt(Parameter.wsize2);
K=max(AssignVec);

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

figure;colormap 'jet'; 
subplot(1,2,1); imagesc(CoOcN);
ylabel ('Co-Occurence matrix');
xlabel(['mean Entropy per row: ',num2str(H)],'Color','red')

subplot(1,2,2);
modes=sum(CoOcN>0,2);
plot(modes);title (strcat('using Thr: ',num2str(Parameter.Spatil.CoOcThr)));
grid on; xlabel('Labels'); ylabel('modes for every cluster')

end

% glcm = graycomatrix(AssignImg,'GrayLimits',[1,K],'NumLevels',K,'Offset',[0,1;1,0],'Symmetric',true);
%CoOc=sum(glcm,3);