function Output=removenoise(Image,Noise,AssignVec)
%% Output=removenoise(Image,Noise,AssignVec)
% 
% Output is the cleaned image, based on the lables in Assign Vec
% Data=Image+Noise (all matricies size row by col )
% if noise is unknown than Image can be Patches , Noise=[]


global Parameter Analysis
AssignVec=AssignVec(:)';
wsize=Parameter.wsize2^0.5;
row=Parameter.row;      col=Parameter.col;
%% BeBug
if isfield(Analysis,'samp') &&  Analysis.DebuggerMode
h(1).fig=figure('Name',['pixel ',num2str(Analysis.samp(1))]);   h(1).i=1;
h(2).fig=figure('Name',['pixel ',num2str(Analysis.samp(2))]);   h(2).i=1;
h(3).fig=figure('Name',['pixel ',num2str(Analysis.samp(3))]);   h(3).i=1;
end
%%
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
temp=zeros(size(X));
W=zeros(1,length(X));

K=max(AssignVec);
%% ##############################################
if Parameter.DeNoise.soft
    disp('Soft assignment Denoising applied')
    l=0;            knn=5;
    Centers=calcCenters(X,AssignVec,K,wsize^2);
    E = X2C (X,Centers,AssignVec,knn);
else
    E = zeros(K,length(X));
    E(sub2ind(size(E),AssignVec,1:length(AssignVec)) )=1;
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

    end
%% ##############################################            
    cluster_0=X(:,AssignVec==k)-center(:,ones(1,clstsize));
    [U, S, V]= svd(cluster_0,'econ');
    SV=S*V';
    UB= U' * B_0;
    SV_hat = shrinkage(SV, 1.1, Parameter.sigma, wsize, 'NULL_LAST_COEFF'); %nulling last coeffs
    UB_hat = shrinkage(UB, 1.1, Parameter.sigma, wsize, 'NULL_LAST_COEFF');
    
    EX(:,(AssignVec==k))= EX(:,(AssignVec==k)) + ones(wsize^2,1)*E(k,AssignVec==k) .*(center(:,ones(1,clstsize))+U*SV_hat);
                                                
    if Parameter.DeNoise.soft
        W(AssignVec==k) = W(AssignVec==k) + E(k,AssignVec==k);
        EX(:,Ind)= EX(:,Ind) + ones(wsize^2,1)*E(k,Ind) .* ( center(:,ones(1,sum(Ind)))+ (U*UB_hat) ); %E(double(k)*ones(wsize^2,1),Ind)
        W(Ind) = W(Ind) + E(k,Ind);
%% BeBug        
        if sum(E(k,Analysis.samp)>0 )&& Analysis.DebuggerMode  %Debug
            temp(:,Ind)= ones(wsize^2,1)*E(k,Ind) .* ( center(:,ones(1,sum(Ind)))+ (U*UB_hat) );
            temp(:,(AssignVec==k))=  ones(wsize^2,1)*E(k,AssignVec==k) .*(center(:,ones(1,clstsize))+U*SV_hat);
            
            f=find(E(k,Analysis.samp)>0);
            for f=f
            pxl=Analysis.samp(f);
            fprintf('pixels %i has a chance of %f to be in cluster %i \n',pxl,E(k,pxl),k)    
            figure(h(f).fig);
            subplot(3,5,5+h(f).i);  imshow(reshape(center,[11,11]),[]);xlabel(sprintf('label %i',k));
            subplot(3,5,10+h(f).i); imshow(reshape(temp(:,pxl),[11,11]),[]);xlabel(sprintf('Pr.= %1.5f',E(k,pxl)));
            h(f).i=h(f).i+1;
            end
        end
%%        
    end
    
    end
end
% fprintf ('the total sum of ewights is %f out of %i pixels\n',sum(W),length(X))        %DEBUG
if  Analysis.DebuggerMode  %Debug
    for f=1:3
        figure(h(f).fig);
        subplot(3,5,2);imshow(reshape(X(:,Analysis.samp(f)),[11,11]),[] );title('Noisy Input');
        subplot(3,5,5);imshow(reshape(EX(:,Analysis.samp(f)),[11,11]),[] );title('Soft Reconsturct');
        subplot(3,5,6);ylabel ('Centers');subplot(3,5,11);ylabel ('Reconstruct')
        if ismatrix(Noise)
            M=floor(wsize/2);
            [j,i]=ind2sub([row-2*M,col-2*M],Analysis.samp(f));
            subplot(3,5,1);imshow(Image(i:i+2*M,j:j+2*M),[] );title('Clean');
        end
    end
end
%% End Debug
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