function ShowCoOc(AssignVec)
global Parameter
row=Parameter.row;      col=Parameter.col;
wsize=sqrt(Parameter.wsize2);
K=max(AssignVec);

AssignImg=col2im(AssignVec,[wsize,wsize],[row,col],'sliding');
[m,n]=size(AssignImg);

glcm = graycomatrix(AssignImg,'GrayLimits',[1,K],'NumLevels',K,'Offset',[0,1;1,0],'Symmetric',true);

CoOc=sum(glcm,3);
CoOcNorm=sum(CoOc,2);
CoOcN=CoOc./CoOcNorm(:,ones(1,size(glcm,1)),:);

H=entropy(CoOc);    HN=entropy(CoOcN);

figure;colormap 'jet'; 
subplot(1,2,1); imagesc(CoOc);axis image
ylabel ('Co-Occurence matrix');
xlabel(['image Entropy: ',num2str(H)],'Color','red')

subplot(1,2,2);imagesc(CoOcN); colorbar;axis image
xlabel(['image Normalize Entropy: ',num2str(HN)],'Color','blue')

end