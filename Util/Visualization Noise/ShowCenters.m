function []= ShowCenters(Centers)
% function []= ShowCenters(Centers)
% shows a sample of the Centres as an image
% uses scaling to be able to see also small patches.

[dim,D,K]=size(Centers);
samp= 25; rows=5;   samp=min([samp,max([K,D])]);%csae less then 25 centers
if K==1;Centers=reshape(Centers,dim,1,K);end


Centers=Centers(:,1,1:samp);
Centers=reshape(Centers,[sqrt(dim),sqrt(dim),samp]);

col=ceil(samp/rows);    %how much columns needed
figure;colormap gray 
for c=1:col             %how much in every row
    m=0;
    while m*col+c<=samp && m<rows
        subplot(rows,col,m*col+c);                                      %[m,c]
        imagesc( imresize(Centers(:,:,m*col+c),[15,15],'bicubic') ,[0,255]);
        axis image off
        m=m+1;
    end

    
end
end

%     subplot(rows,col,c);                                            %[1,c]
%     imagesc( imresize(Centers(:,:,c),[15,15],'bicubic') );
%     subplot(rows,col,col+c);                                        %[2,c]
%     imagesc( imresize(Centers(:,:,col+c),[15,15],'bicubic') );
%     subplot(rows,col,2*col+c);                                      %[3,c]
%     imagesc( imresize(Centers(:,:,2*col+c),[15,15],'bicubic') );
%     subplot(rows,col,3*col+c);                                      %[4,c]
%     imagesc( imresize(Centers(:,:,3*col+c),[15,15],'bicubic') );
%     subplot(rows,col,4*col+c);                                      %[5,c]
%     imagesc( imresize(Centers(:,:,4*col+c),[15,15],'bicubic') );