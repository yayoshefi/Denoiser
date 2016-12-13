function Output=removenoise(Image,Noise,AssignVec)
%% Output=removenoise(Image,Noise,AssignVec)
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
if length (AssignVec)~=length(X)  %FIX: cause of window on Label, what image 2 reconstruct with 0 border?
    mprime=sqrt(length(AssignVec)); m=row-wsize+1;
    %mprime=Analysis.LabelSize(1)
    padding=floor( (m-mprime)/2 );
    AssignVec=padarray(reshape(AssignVec,[mprime,mprime]),[padding,padding]);
    AssignVec=AssignVec(:)';
    disp ('something went wrong, AssignVec was not full size');
end 
EX=zeros(size(X));
W=zeros(1,length(X));

K=max(AssignVec);
%% ##############################################
if Parameter.DeNoise.soft
    disp('Soft assignment Denoising applied')
    l=0;            knn=5;
    Centers=calcCenters(X,AssignVec,K,wsize^2);
    E = X2C (X,Centers,AssignVec,knn);
end
%%  ##############################################
B_0=double.empty(wsize^2,0);

for k=1:K
    clstsize=sum(AssignVec==k); clearvars unsort
    if Parameter.ORACLE
            center=Analysis.ORACLE.Centers(:,:,k);
    else    center=mean(X(:,AssignVec==k),2);
    end
%% robust is a way to handle outliers in the decompostion of the cluster
    if Parameter.DeNoise.Robust
        if k==1; disp('Robust Denoised applied');end
        cluster=X(:,AssignVec==k);
        pointdist = sum(cluster.^2) - 2*center'*cluster + center'*center;
        clusterDist = abs(pointdist).^0.5;
        [~,I]=sort(clusterDist);
        unsort(I)=1:clstsize;
        
        New_size= round(clstsize*0.92);
        if New_size==0; New_size=1; end
%         if New_size<Parameter.wsize2        %DeBugging
%             fprintf('\nlabel %i the New cluster holds %i data points',k,New_size)
%         end
        A = cluster( :,I(1:New_size) );         % A - is the new (closer) cluster
        B = cluster( :,I(New_size+1:end) );     % B - is the outliers from A
        
        New_center   = mean(A,2);
        A_0 = A - New_center(:,ones(1,New_size));
        B_0 = B - New_center(:,ones(1,clstsize-New_size));

        [U, S, V]= svd(A_0,'econ');
        SV= S*V';           %A_hat = A-SV;
        UB= U'*B_0;         %B_hat = B-UB;
        SV_hat=shrinkage([SV,UB], 1.1, Parameter.sigma, wsize, 'NULL_LAST_COEFF'); %nulling last coeffs
        EX(:,AssignVec==k)=New_center(:,ones(1,clstsize))+U*SV_hat(:,unsort);
    else
%% Normal decomposition, using all patches
    
%% ##############################################        
    if Parameter.DeNoise.soft
        if k==1;disp('Soft assignment Denoising applied');end
        Ind = ( (E(k,:) > 0) & ~(AssignVec==k) );
        B   = X(:,Ind );
        B_0 = B-center(:,ones(1,size(B,2)));
%         i want to put this box inside the regular part, i need to update:
%         setting the choosen data points to all those where ( E[k,:] > 0 )
%         and using an aggregate value to EX.    EX(:,ind)= EX(:,ind)+.... E[k*ones(length(center),1),:]* Recosntruct
    end
%% ##############################################            
    cluster_0=X(:,AssignVec==k)-center(:,ones(1,clstsize));
    [U, S, V]= svd(cluster_0,'econ');
    SV=S*V';
    UB= U' * B_0;
    SV_hat = shrinkage(SV, 1.1, Parameter.sigma, wsize, 'NULL_LAST_COEFF'); %nulling last coeffs
    UB_hat = shrinkage(UB, 1.1, Parameter.sigma, wsize, 'NULL_LAST_COEFF');
    
    EX(:,(AssignVec==k))= EX(:,(AssignVec==k)) + center(:,ones(1,clstsize))+U*SV_hat;
    
    if Parameter.DeNoise.soft
        W(AssignVec==k) = W(AssignVec==k) + E(k,AssignVec==k);
        EX(:,Ind)= EX(:,Ind) + ones(wsize^2,1)*E(k,Ind) .* ( center(:,ones(1,sum(Ind)))+ (U*UB_hat) ); %E(double(k)*ones(wsize^2,1),Ind)
        W(Ind) = W(Ind) + E(k,Ind);
    end
    
    end
end
fprintf ('the total sum of ewights is %f out of %i pixels',sum(W),length(X))        %DEBUG

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


function [Centers] = calcCenters(Data,AssignVec,K,wsize2)
Centers=zeros(wsize2,1,K);
for k=1:K
   Centers(:,1,k)=mean ( Data(:,AssignVec==k),2 ) ;
end
end

function [E,AssignMat] = X2C (Data,Centers,AssignVec,knn)
[wsize2,pnum]=size (Data);
[K]=size(Centers,3);
E=zeros(K,pnum);
sigma=wsize2/4;

[AssignMat, Distances]=Dist2SubSpace (Data,Centers,'knn',knn);
%% Dist2Subspace
% AssignMat=zeros(knn,pnum);
% Distances=inf*ones(knn,pnum);
% D=inf*ones(K,pnum);
% 
% h=waitbar(0,'Computing data distances..','CreateCancelBtn','setappdata(gcbf,''cancel'',1)');
% setappdata(h,'cancel',0);
% for k=1:K
%     if getappdata(h,'cancel')
%     break
%     end
%     waitbar(k/K);
%     pointdist=sum(Data.^2)-2*Centers(:,:,k)'*Data+Centers(:,:,k)'*Centers(:,:,k);
%     if Centers(1,:,k)==inf; pointdist=inf; end %fix empty clusters
%    tmpdist=(abs(pointdist)).^0.5;
% 
%    if knn==1
%        track=tmpdist<Distances(end,:);
%        Distances(end,track)=tmpdist(track);
%        AssignMat(end,track)=k;
%    else
%        D(k,:)=tmpdist;
%    end
% end
% if knn>1
%    waitbar(k/(K+10),h,'Sorting smallest distances ...');
%    [Distances,I]=sort(D);
%    Distances=Distances(1:knn,:);
% 
%    %sort the matrix AssignVec by distances
%    [col_sub,AssignMat]=meshgrid(1:pnum,1:K);
%    sort_indx=sub2ind([K,pnum],I(:),col_sub(:));
%    AssignMat=reshape(AssignMat(sort_indx),K,[]);
%    AssignMat=AssignMat(1:knn,:);
% end
% delete(h);
%%

col_sub=meshgrid(1:pnum,1:knn);
ind=sub2ind([K,pnum],AssignMat(:),col_sub(:));

E(ind)=exp(Distances(:)/(-sigma^2));

cumE=sum(E,1);
NormalizeFactor=1./cumE';
NormalizeFactor( cumE'==0 )=0;                  %outlier point- far from all centers
if sum (cumE==0)>log(pnum)
    error ('AffinityFunc:Outliers','%g points where Sum[Pr(k)]=0 out of %g total points',sum (cumE==0),pnum)
end
E=DiagonalMult(E,1.*NormalizeFactor,'r');           %normalize the affinity matrix

end