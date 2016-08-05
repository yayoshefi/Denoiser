%% function [CoOc,P_Ni] =  lcm(arg1,AssignMethod,Type)
% calculates the Co-Occurrence matrix
% Input:

% AssignMethod - can be 'hard' or 'soft' assignment
% arg1 -    case using soft arg1= AssignVec or AssignImg
%           case using hard arg1= Pi (3rd dim specifing the affinity
%                                      to each center) 
% Type choosed between 
%           'CC' - conditional CC
%           'M'  - Mutual information 
%           'JP' - Joint Probabilty


function [CoOc,P_Ni]=lcm(arg1,AssignMethod,Type)
if ~exist('Type','var');Type='CC';
elseif isempty(Type); Type='CC'; end

global Parameter Analysis
wsize=sqrt(Parameter.wsize2);
m=Parameter.row-wsize+1;    n=Parameter.col-wsize+1;
M=m*n;
NN=Parameter.spatial.NN;  %window2
padding=floor(NN/2);
K=Analysis.K;


filterType='gaussian';      %  'gaussian' or 'rect'
[NormailzeFactor,filt]= SetNeigborWindow(NN,m,n,filterType);   

switch lower(AssignMethod)
    case 'hard'
        if size (arg1,3)>1
            error('incosiset size for assignVec, Third dim is %u',size (arg1,3));
        elseif size(arg1,1)==1|| size(arg1,2)==1
             AssignImg=col2im(arg1,[wsize,wsize],[Parameter.row,Parameter.col]);
        else AssignImg=arg1;    arg1=AssignImg(:);
        end

        Indicator=sparse(double(arg1(:)),1:M,ones(1,M),K,M);
        switch filterType
            case 'rect'
                AssignImg=padarray(AssignImg,[padding,padding],-1);
                Neigbour=im2col(AssignImg,[NN,NN],'sliding');
                Neigbour(ceil(NN^2/2),:)=[];
                Hstar=histc(Neigbour,1:K,1)';        
                H= bsxfun(@rdivide,Hstar,NormailzeFactor);
            case 'gaussian'
                P_i=reshape(full(Indicator)',[m,n,K]);
                P_Ni=reshape( convn(P_i,filt,'same'),[M,K] );       %H*i
                H= bsxfun(@rdivide,P_Ni,NormailzeFactor);   %Hi
        end
        C= Indicator*H;
        
        clstCnt=histcounts(arg1(:),1:K+1);
        NullClusters=(clstCnt==0);
        Pl=clstCnt/sum(clstCnt);
         
        CC=diag(1./clstCnt)*C;                  CC(NullClusters,:)=0;
        switch Type
            case 'CC'
                CoOc=CC;
            case 'MI'
                CoOc=CC*diag(1./clstCnt);       CoOc(:,NullClusters)=0;
                %MI=  JP/( Pl(l) Pl(k) )
            case 'JP'
                CoOc=C/M;
            case 'PMI'
                CoOc=diag(1./Pl) * (C/M).^2 * diag(1./Pl);
        end
        P_Ni=H;
        
    case 'soft'
        if size (arg1,3)==1
            error('incosiset size for P_i ,input should be 3 dim matrix\n instead the size is: %i by %i', size(arg1,1),size(arg1,2) );
        elseif size(arg1,1)==1|| size(arg1,2)==1
             P_i=reshape(arg1,[m,n,size(arg1,3)]);
        else P_i=arg1;
        end
        P_Ni=reshape( convn(P_i,filt,'same'),[M,K] );       %H*i
        P_Ni_mat= bsxfun(@rdivide,P_Ni,NormailzeFactor);   %Hi
        Pi_mat  = reshape(P_i  ,[M,K])';

        Pl=mean(Pi_mat,2);
        JP=(Pi_mat*P_Ni_mat/M);
        switch Type
        case 'CC'
            CoOc=diag(1./Pl) * JP ;
        case 'MI'
            hl=sum(Pi_mat,2);
            CoOc=diag(1./hl) * JP * diag(1./hl);
        case 'PMI'
            CoOc=diag(1./Pl) * JP.^2 * diag(1./Pl);
        case 'JP'
            CoOc=JP;
        end    
end

end

function [NormailzeFactor,filt]= SetNeigborWindow(NN,m,n,type)
if ~exist('type','var'); type='rect';end
sigma=1;
switch type
    case 'rect'
        filt=ones(NN);
        filt( (NN+1)/2,(NN+1)/2 )=0;
    case 'gaussian'
        filt=fspecial(type,NN,sigma);
        filt( (NN+1)/2,(NN+1)/2 )=0;        
end
NormailzeFactor=conv2(ones(m,n),filt,'same');
NormailzeFactor=NormailzeFactor(:);             %used for bounderies 
end

function [Out]=Convert (Input,convert_from,m,n,K)
% used for debugging when comparing between hard and soft
switch convert_from
    case 'soft'     % soft 2 hard
        Pi=Input;
        AssignImg=zeros (m,n);
        for k=1:K
            AssignImg=AssignImg+k*Pi(:,:,k);
        end
        Out=AssignImg;
    case 'hard'     % hard 2 soft
        AssignImg=Input;
        Pi=single(zeros(m,n,K));
        for k=1:K
        Pi(:,:,k)=(AssignImg==k);
        end
        Out=Pi;
end
end

