function [AssignVec2]=SpatialContext (Data,AssignVec,Centers,varargin)
%% [AssignVec2]=SpatialContext (Data,AssignVec,Centers,varargin)
%
% find the lables for each path using Image context
% using spectral clustring or graph cut (#Shai bagon)

% [AssignVec2]=SpatialContext (...,'s',S,'basis'Basis)
%   as in the rest of function
%  for computing diffrent distance metric- mahalanobis


global Parameter Analysis

% ----- Parameters -----
Parameter.metric ='euclidean';
K=size(Centers,3);
wsize=sqrt(Parameter.wsize2);
m=Analysis.LabelsSize(1);
n=Analysis.LabelsSize(2);


pnames = {'s'           'basis' };
dflts =  {ones(wsize^2,1,K),0 };
[Energy,Basis] = internal.stats.parseArgs(pnames, dflts, varargin{:});

switch Parameter.Context
    case 'spectral'
%         if Parameter.Spatil.lambda==0; Parameter.Spatil.lambda=0.1; end
%         E=CoherentAffinity (Patches, Centers,'s',Energy,'basis',Basis,'lambda',lambda);
        [E]=affinity (Data, Centers,20,1);
        B = Spatialaffinity ();
        E_spatial=[E;Parameter.Spatil.lambda*B];
        [AssignVec2]=MySpectralClustering(E_spatial,'e');
                
    case 'graphcut' %mrf solver 
        [E]=affinity (Data, Centers,0,1);
        Dc=exp(reshape(E',[m,n,K]).*10^2);
        Sc = ones(K) - eye(K);
        
        [Hc, Vc] = SpatialCues( reshape(Data(1,:),[m,n]) );
        
        gch = GraphCut('open', Dc, Parameter.Spatil.lambda*Sc, exp(-Vc*5), exp(-Hc*5) );
        [gch, L] = GraphCut('expand',gch,5);
        gch = GraphCut('close', gch);
        
        AssignVec2=L(:)';
    case 'rl' %relaxation labling
        knn=0; value=true;
        [E]=affinity (Data, Centers,knn,value);          %,varargin{1:6});
        [E_directional,comp]=Compatability (E,AssignVec);
        S=sum( E_directional.*comp(ones(K,1),:,:),3 );
        
        p=logical(padarray(ones(m-2,n-2),[1,1]));
        p=p(:)';
        Enew=E(ones(K,1),p).*(1+S);
        [Pr, AssignVec2]=max(Enew);
        
    case 'mrf'
            [AssignVec2]=MRFContext(Data,AssignVec,Centers);
    case 'entropy'  %minimize Co-Occurence sparsity
        %use delta function to converge to this distriboution
        AssignImg=col2im(AssignVec,[wsize,wsize],[Parameter.row,Parameter.col]);
        NN=Parameter.Spatil.NN;
        Neigbour=im2col(AssignImg,[NN,NN],'sliding');
        Neigbour(ceil(NN^2/2),:)=[];
        LocalHist=histc(Neigbour,1:K,1);        
        
        Histoids=K*eye(K);
        [AssignVec2,C]=kmeans(LocalHist',K,'start',Histoids,'distance','correlation');
        
    case 'comeans'
        [S]=affinity (Data, Centers,0,true)';
        Lhat=AssignVec;
        NN=Parameter.Spatil.NN;  %window2
        padding=floor(NN/2);
        Analysis.LabelsSize=Analysis.LabelsSize-(NN-1);
        m=Analysis.LabelsSize(1)    ;n=Analysis.LabelsSize(2);
        p=logical(padarray(ones(m,n),[padding,padding]));
        
        Analysis.LabelsSize=Analysis.LabelsSize+(NN-1); %restore values to origin
        
        samp=randperm(size(S,1),3);
        for iter=1:4  
            AssignImg=col2im(Lhat,[wsize,wsize],[Parameter.row,Parameter.col]);
%             AssignImg=padarray(AssignImg,[padding,padding],-1);
            
            Neigbour=im2col(AssignImg,[NN,NN],'sliding');
            Neigbour(ceil(NN^2/2),:)=[];
            H=histc(Neigbour,1:K,1)';
            
            Indicator=sparse(Lhat(p),1:m*n,ones(1,m*n),K,m*n);
            CC=Indicator*H;
            CCNorm=sum(CC,2);
            CCN=CC./CCNorm(:,ones(1,K),:); CCN(CC==0)=0; %to avoid 0/0=nan
            
            Parameter.Spatil.CoOcThr=0.01;
            
            CCN(CCN<Parameter.Spatil.CoOcThr)=0;
            CCNorm=sum(CCN,2);
            CCthr=CCN./CCNorm(:,ones(1,K),:); CCthr(CCN==0)=0;
            
            L=S;
            L(p,:)=(1-Parameter.Spatil.lambda)*reshape(S(p,:),m*n,K)+Parameter.Spatil.lambda/(NN^2-1)*H*CCthr;
            [Pr,Lhat]=max(L,[],2);
            
            Debug (CCthr,Lhat,Pr,iter);
            ShowProb (L,samp);

            
            [Centers,~,Lhat,~,~]=UpdateCenter(Data,Lhat,false);
            [S]=affinity (Data, Centers,0,true)';
        end
        
        AssignVec2=Lhat;
end
 if length(AssignVec2)~=length(AssignVec) % so that AssignVec2 will be always same size
     mprime=sqrt(length(AssignVec2)); m=Parameter.row-wsize+1;
     %mprime=Analysis.LabelsSize(1)
     padding=floor( (m-mprime)/2 );
     p=logical(padarray(ones(sqrt(length(AssignVec2))),[padding,padding]));
     AssignVec(p)=AssignVec2;
     AssignVec2=AssignVec;

 end
end

function [hC vC] = SpatialCues(im)
g = fspecial('gauss', [13 13], sqrt(13));
dy = fspecial('sobel');
vf = conv2(g, dy, 'valid');
sz = size(im);

vC = zeros(sz(1:2));
hC = vC;

for b=1:size(im,3)
    vC = max(vC, abs(imfilter(im(:,:,b), vf, 'symmetric')));
    hC = max(hC, abs(imfilter(im(:,:,b), vf', 'symmetric')));
end
end

function B = Spatialaffinity ()
global Parameter

wsize=sqrt(Parameter.wsize2);
m=Parameter.row-wsize+1;
n=Parameter.col-wsize+1;
[ColI,RowI]=meshgrid(1:n,1:m);  % coordinates of an the image

switch Parameter.Spatil.spatialdist
    case 'landmarks'
        L=floor((m*n)^0.35);     %number of land marks
        sigma=500;

        Row=1:Parameter.row;
        Col=1:Parameter.col;

        LandMarks=[Parameter.row*rand(1,L) ; Parameter.col*rand(1,L)];
        signiture=zeros(m*n,L);
        for l=1:L
            signiture(:,l)=sum([ColI(:),RowI(:)].^2,2) - 2*[ColI(:),RowI(:)]*LandMarks(:,l) + ...
                LandMarks(:,l)'*LandMarks(:,l);
        end
        signiture=sqrt(signiture'); % 2 norm
        B=exp(-signiture/sigma);
        Bn=sum(B);
        B=B./Bn(ones(L,1),:);
        
    case 'simplenoramlize'
        B=[RowI(:)';ColI(:)'];
        Bn=sum(B);
        B=B./Bn([1;1],:);
        
    case 'decomposition'
        sigma=Parameter.Spatil.sigma;
        l=1:max(m,n);
        l2=repmat(l.*l,length(l),1);  %n^2
        D2= l2+l2'-2*l'*l;  % D(i,j)= (n(i)-n(j))^2
        D2(D2>sigma^2)=sigma^2;         %Tukey trancation
        gxy=exp(-D2./sigma^2);  % Gaussian function


        [U,S,V]=svd(gxy);  %decompose 1D Gaussian Function g(|xi-xj|)

        SS=cumsum(diag(S)); %cumulative sum
        SS=SS./max(SS);
        I=find(SS>0.90,1);  %finds the number of eigenvectors spanning 90% 
        B1D=U(:,1:I)*sqrt(S(1:I,1:I));
%         agxy=B*B';
        BB_col=B1D(ColI(:),:);
        BB_row=B1D(RowI(:),:);
        
        % multiply all eigenvectors pairs
        L=size(BB_col,2);
        B=[];
        for l=1:L
            BBk=repmat(BB_col(:,l)',L,1);
            B=[B; BBk.*BB_row'];
            
        end;
        Bn=sum (B.^2,2);
        B=DiagonalMult(B,Bn,'l');       % normalize proximity        
    case 'none'
        B=[];        
end
Parameter.B=B;        %De-bugging

end

function[] = Debug (CoOc,Lhat,Pr,iter)
global Parameter
wsize=sqrt(Parameter.wsize2);

% ShowCoOc(Lhat); set(gcf,'Name',strcat('Co-Occurence ',num2str(iter), ' iteration'));
Pr_Img=col2im(Pr,[wsize,wsize],[Parameter.row,Parameter.col]);
tmp_Labels=col2im(Lhat,[wsize,wsize],[Parameter.row,Parameter.col]);

figure('Name',strcat('temporal image properties ',num2str(iter), 'iteration'));
subplot(2,2,3);
imagesc(CoOc);colormap jet; title ('Co-Occurence matrix');

LCoOc=log2(CoOc);
LCoOc(CoOc==0)=0;  % no nan resulting from inf*0;
H_row=-sum( CoOc.*LCoOc,2 );
H=mean(H_row);

xlabel(strcat('mean entropy for each row is: ',num2str (H)));
subplot(2,2,4);
modes=sum(CoOc>0,2);
plot(modes);title (strcat('using Thr: ',num2str(Parameter.Spatil.CoOcThr)));
axis([1,size(CoOc,1),0,15]);grid on
xlabel('Labels');ylabel('modes for every cluster')

subplot(2,2,1);
imagesc(Pr_Img);title ('Probabilty to be in Lhat'); colormap jet
xlabel(strcat('mean prob.= ',num2str(mean(Pr)) ))
subplot(2,2,2)
imagesc(tmp_Labels);title ('Temp Label image'); colormap jet
xlabel(strcat(num2str(length(unique(Lhat))) ,' diffrent labels'))

end