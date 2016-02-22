function [AssignVec]=MySpectralClustering(Data,affinity,varargin)
%% MySpectalClustering(Data,K,affinity)
%
% Spectral Clustring is based on the similarity between signal
% the Similarity matrix is built from the Affinity matrix E
% Sim=E*E'    E is [pnum by d] matrix
% 
%
% ---- Input ---
% K is the number of clusters to be found.
%
% affinity chooses how does the affinity matrix is obtained
% affinity can be: 'ksvd' or 'e' or 'nystrom' or 'dist'
%
% 'ksvd'- is for bulding the dictionary,
% using the sparse represantation as the affinity matrix.
% the function uses sparsity and dictsize as params
% E = Alpha'
%
% 'e'-is if Data is allready an affinity matrix such that W=E*E'.
% Data need to be [d by pnum]
%
% 'nystrom' - method is based on LandMarks and rest of the points
% Sim=[A  B
%      B' C]
%such that A is the full Similarity given for small data points (say centers)
% and B is the Affinity between the patchs and Centers.
%
% [AssignVec]=MySpectralClustering(...,'Centers',Centers)
% Centers can be a matrix [d X K] or 3d array [d X 1 X K] 
% K is the number of LandMark, K<< pnum
%
% 'dist'- is when the between data points distance is given.
% the Similarity matrix is built using gaussian func to the dist
% ###### currently don't work because of memory problems ##########

 
tic 
% Parameters
global Parameter
K=Parameter.Spectral.clustrsNUM;
pnames = {'Centers'};
dflts =  {[]};
[Centers]= internal.stats.parseArgs(pnames, dflts, varargin{:});
pnum=size(Data,2);


%% initializtion

switch affinity
    case 'ksvd'
        signum=size(Data,2);
        samp=randperm(signum,round(Parameter.subsample*signum));
        SampData=Data(:,samp);

        PARAMS=struct('data',SampData,'Tdata',Parameter.Spectral.sparsity,...
            'dictsize',Parameter.Spectral.dictsize,'iternum',4);

        [Dict,~]=ksvd(PARAMS,'');
        disp(['end of ksvd: ',num2str(round(toc)),' sec'])
        Alpha=omp(Dict,Data,Dict'*Dict,Parameter.Spectral.sparsity);
        clearvars Dict
        
        if Parameter.Spectral.HardThr
             fh=@logical;
        else fh=@abs;       end  
        
        Alpha=spfun(fh ,Alpha);
        D=calcD(Alpha');
        D1_2=D.^(-0.5);
        B=sparse(DiagonalMult(Alpha,D1_2','r'));
        
    case 'e'
%         E=Data';
        D=calcD(Data');
        D1_2=D.^(-0.5);           %D^-(0.5)
        B=DiagonalMult(Data,D1_2','r');

        
    case 'dist'
        % compute the similarity matrix, maybe by using correlation using
        % full distances between Point2Center
        % or Sim= Similarity2 (Data,10,1);
    case 'nystrom'
        %W=[A  B
        %   B' C]
        % B==Data is all the Centers affinity to LM , B=E [K by pnum]
        
        Centers=squeeze(Centers);
        K=size(Centers,2);
        A=zeros(K);
        for k=1:size(Centers,2)
            dist=sum(Centers.^2)-2*Centers(:,k)'*Centers+Centers(:,k)'*Centers(:,k);
            A(k,:)=sqrt(abs(dist));
        end
        % normalization - DIDNT UNDERSTAND
        d1=sum([A;Data'],1);
        d2=sum(Data,1)+sum(Data',1)*pinv(A)*Data;
        dhat=sqrt(1./[d1, d2]');
        A=A.*(dhat(1:K)*dhat(1:K)');
        Data=Data.*(dhat(1:K)*dhat(K+(1:pnum))');
%         [U,G,~]=svd(A,'econ');
        Asi=sqrtm(pinv(A));
        S=A+Asi*(Data*Data')*Asi;
        [Us,Gs,~]=svd(S,'econ');
        
% ******** nystrom method needs large enough A ************
% V= [A;Data']*Asi*Us*( Gs^(-1/2) )
% B= sqrt(Gs) * V'
% W_hat=B'*B
        B=(sqrt(Gs))*Us'*Asi'*[A];      %B in respect to 1st part==A
        band=1000;
        for bp=1:band:pnum
            chunk=bp:min(bp+band-1,pnum);
            B=[B,(sqrt(Gs))*Us'*Asi'*[Data(:,chunk)] ];
        end
        
end
% H=D1_2*Sim*D^-(0.5);
% B=E'*D^-0.5
% B is [K by n]  K<<n

K=min(K,size(B,1));
if strcmp(affinity,'ksvd') ;[U,S,~]=svds(B*B',size(B,1));
else [U,S,~]=svd(B*B'); end
disp(['end of svd decomposition: ',num2str(round(toc)),' sec'])

S=PseudoDiagInv(diag(S));
V=DiagonalMult(B'*U,(S.^0.5),'r');

X=V(:,1:K);
Xnorm=( sum(X.^2,2) ).^0.5;
Xn=DiagonalMult(X,1./Xnorm,'l');
% Xn=X./Xnorm(:,ones(1,K));
disp(['before kmeans: ',num2str(round(toc)),' sec'])

if strcmp(Parameter.Spatil.spatialdist,'none') && strcmp(Parameter.Context,'spectral')  
    % case we use spatial context and haven't inserterd none in the
    % affinity matrix
    m=Parameter.row-sqrt(Parameter.wsize2)+1;
    n=Parameter.col-sqrt(Parameter.wsize2)+1;
    Parameter.Spatil.lambda=1e-2;
    disp(['using SLIC for context clustring. Spatial Lambda: ',num2str(Parameter.Spatil.lambda)])
    SLIC=1;%sqrt(m*n/K);
    [ColI,RowI]=meshgrid(1:n,1:m);  % coordinates of an the image
    Xn=[Xn,Parameter.Spatil.lambda*SLIC*RowI(:),Parameter.Spatil.lambda*ColI(:)];
end

opts = statset('Display','off','MaxIter',50);
if Parameter.Spectral.Fast
    samp=randperm(pnum,round(Parameter.subsample*1.4*pnum));
    SampXn=Xn(samp,:);
    [~,C]=kmeans(SampXn,K,'onlinephase','off','Start','cluster','options',opts);
    PseudoCenters=reshape(C',[size(Xn,2),1,K]);
    [AssignVec,~]=Dist2SubSpace(Xn',PseudoCenters);
else
    AssignVec=kmeans(Xn,K,'onlinephase','off','Start','cluster','options',opts);%,'Replicates',2
end
disp(['end of k means: ',num2str(round(toc)),' sec'])
end


%{
1) E        -affinity matrix [E=Alpha' or Data']
2) W=E*E'    -weights matrix
    L=D-W;  H=D^-0.5*W*D^-0.5;  L(sym)=I-H       L-Laplcian
H=B'*B;         [B=E'*D^-0.5]

[U,S,V]=svd(B)
H=V*S^2*V
V=?

[U,S^2,U]=svd(B*B')
V=B'*U*S^-1

Choose k largest eigenvectors of V
normalize each row to be ||Vi||=1
kmeans row-wise
each point is in the same cluster as the matching row is V
%}

