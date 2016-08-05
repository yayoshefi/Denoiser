function []=BgFg (AssignImg,CoOc)

[r,c]=size(AssignImg);
K=size(CoOc,1);

Like=zeros(r-2,c-2);

R=AssignImg(:,2:end);
L=AssignImg(:,1:end-1);

U=AssignImg(1:end-1,:);
D=AssignImg(2:end,:);
%Right
LikeR=reshape(CoOc(sub2ind([K,K],R(:),L(:))),size(R));
Like=Like+LikeR(2:end-1,1:end-1);
%Left
LikeL=reshape(CoOc(sub2ind([K,K],L(:),R(:))),size(L));
Like=Like+LikeL(2:end-1,2:end);

%Up
LikeU=reshape(CoOc(sub2ind([K,K],U(:),D(:))),size(U));
Like=Like+LikeU(2:end,2:end-1);
%Down
LikeD=reshape(CoOc(sub2ind([K,K],D(:),U(:))),size(D));
Like=Like+LikeD(1:end-1,2:end-1);


figure;
imshow((Like).^2,[]);
end