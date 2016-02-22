function []= ShowCenters(Centers)
% function []= ShowCenters(Centers)
% shows a sample of the Centres as an image
% uses scaling to be able to see also small patches.

[dim,~,K]=size(Centers);
samp= 25; rows=5;
if K==1;Centers=reshape(Centers,dim,1,K);end


Centers=Centers(:,1,1:samp);
Centers=reshape(Centers,[sqrt(dim),sqrt(dim),samp]);

col=ceil(samp/rows);
figure;colormap gray 
for c=1:col
    subplot(rows,col,c);
    imagesc( imresize(Centers(:,:,c),[15,15],'bicubic') );
    subplot(rows,col,col+c);
    imagesc( imresize(Centers(:,:,col+c),[15,15],'bicubic') );
    subplot(rows,col,2*col+c);
    imagesc( imresize(Centers(:,:,2*col+c),[15,15],'bicubic') );
    subplot(rows,col,3*col+c);
    imagesc( imresize(Centers(:,:,3*col+c),[15,15],'bicubic') );
    subplot(rows,col,4*col+c);
    imagesc( imresize(Centers(:,:,4*col+c),[15,15],'bicubic') );
end
end