function [Output]=NonMaxSus(Data,Centers,Input)
global Parameter
row=Parameter.row;  col=Parameter.col;
wsize=sqrt(Parameter.wsize2);
K=Parameter.values.kmeans;
SuppressionWindow=Parameter.spatial.NN;
SuppressionWindow=5;

[AssignVec, Distances]=Dist2SubSpace (Data,Centers);

AssignImg=col2im(AssignVec,[wsize,wsize],[row,col]);
DistImg=col2im(Distances,[wsize,wsize],[row,col]);


MinDist=ordfilt2(DistImg,1,true(SuppressionWindow),'symmetric');

SuppressedPoints= (MinDist==DistImg);

SuppressedLabels=AssignVec(SuppressedPoints(:));

SuppressedData=Data(:,SuppressedPoints(:)');

[SuppAssignVec,C]=kmeans(SuppressedData',K,'onlinephase','off');
SuppCenters=reshape(C',[size(Data,1),1,K]);

% [AssignVec,~]=Dist2SubSpace(Data,Centers,'dim',zeros(1,1,Parameter.values.kmeans));

%% REMOVE NOISE
if exist('Input','var')
    X=SuppressedData;
    EX=zeros(size(X));

    for k=1:K
        clstsize=sum(SuppAssignVec==k);
    %     center=SuppCenters;
        center=mean(X(:,SuppAssignVec==k),2);
        cluster_0=X(:,SuppAssignVec==k)-center(:,ones(1,clstsize));
        [U, S, V]= svd(cluster_0,'econ');
        SV=S*V';

        SV=shrinkage(SV, 1.1, Parameter.sigma, wsize, 'NULL_LAST_COEFF'); %nulling last coeffs

        EX(:,SuppAssignVec==k)=center(:,ones(1,clstsize))+U*SV;
    end
    window=ones(wsize);
    %window=fspecial('gaussian',wsize2^0.5);

    Output=zeros(row,col);
    Weights=zeros(row,col);

    index=find(SuppressedPoints);
    for p=1:length(index)
        patch=reshape(EX(:,p),wsize,wsize);
        [i,j]=ind2sub([row-wsize+1,col-wsize+1],index(p));
        Output(i:i+wsize-1,j:j+wsize-1)=Output(i:i+wsize-1,j:j+wsize-1)+patch.*window;
        Weights(i:i+wsize-1,j:j+wsize-1)=Weights(i:i+wsize-1,j:j+wsize-1)+window;
    end

    Output=Output./Weights;

    % fill in gaps
    gaps=isnan(Output);
    Output(gaps)=Input(gaps);

else % no reconstruction only the new centers
    Output = SuppCenters;
end

end