function Output=removenoise(Image,Noise,AssignVec)
%% Output=removenoise(Image,Noise,AssignVec)
%
% OLD VER- No protection in case Assignvec2 < AssignVec
%
% Output is the cleaned image, based on the lables in Assign Vec
% Data=Image+Noise (all matricies size row by col )
% if noise is unknown than Image can be Patches , Noise=[]


global Parameter Analysis

wsize=Parameter.wsize2^0.5;
row=Parameter.row;      col=Parameter.col;

if ismatrix(Noise)
    
    X=im2col(double(Image)+Noise,[wsize,wsize],'sliding');
else
    X=Image;
end

EX=zeros(size(X));

K=max(AssignVec);

for k=1:K
    clstsize=sum(AssignVec==k);
    center=mean(X(:,AssignVec==k),2);
    cluster_0=X(:,AssignVec==k)-center(:,ones(1,clstsize));
    cluster=X(:,AssignVec==k);
    [U, S, V]= svd(cluster_0,'econ');
    SV=S*V';
    
    SV=shrinkage(SV, 1.1, Parameter.sigma, wsize, 'NULL_LAST_COEFF'); %nulling last coeffs
    
    EX(:,AssignVec==k)=center(:,ones(1,clstsize))+U*SV;
end
if ismatrix(Noise) && Analysis.Show
    P=im2col(Image,[wsize,wsize],'sliding');
    Analysis.PatchRmse=sqrt( mean((P-EX).^2) );
    Analysis.DC_Change=mean(X)-mean(EX);     %DC of noisy signal - DC after Clean
    
%     ShowCleaning (Image,Noise,EX);
end

window=ones(wsize);
%window=fspecial('gaussian',wsize2^0.5);

Output=zeros(row,col);
Weights=zeros(row,col);

for p=1:length(X)
    patch=reshape(EX(:,p),wsize,wsize);
    [i,j]=ind2sub([row-wsize+1,col-wsize+1],p);
    Output(i:i+wsize-1,j:j+wsize-1)=Output(i:i+wsize-1,j:j+wsize-1)+patch.*window;
    Weights(i:i+wsize-1,j:j+wsize-1)=Weights(i:i+wsize-1,j:j+wsize-1)+window;
end

Output=Output./Weights;

end