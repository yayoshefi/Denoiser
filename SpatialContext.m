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
wsize=sqrt(Parameter.wsize2);
if ~isfield(Analysis,'samp');   Analysis.samp=randperm(length(AssignVec),3);    end

switch Parameter.Context
    case 'spectral'
        AssignVec2=SpectralContext (Data,AssignVec,Centers,varargin);                
    case 'graphcut'
        AssignVec2=MRFContext_addonSolver(Data,AssignVec,Centers);
    case 'rl'
        AssignVec2=RelaxLabel(Data,AssignVec,Centers);
    case 'mrf'
        [AssignVec2]=MRFContext(Data,AssignVec,Centers);
    case 'entropy'  
        AssignVec2= MinEntropyCC (Data,AssignVec,Centers);
    case 'mutualdist'
        AssignVec2=MutualDist(Data,AssignVec,Centers);
    case 'comeans'
         AssignVec2=CoMeans(Data,AssignVec,Centers);
end
 if length(AssignVec2)~=length(AssignVec) % so that AssignVec2 will be always same size
     mprime=sqrt(length(AssignVec2)); m=Parameter.row-wsize+1;
     %mprime=Analysis.LabelsSize(1)
     padding=floor( (m-mprime)/2 );
     p=logical(padarray(ones(sqrt(length(AssignVec2))),[padding,padding]));
     AssignVec(p)=AssignVec2;
     AssignVec2=AssignVec;

 end
 Analysis.K2=length(unique(AssignVec2));
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
imagesc(log(CoOc+1));colormap jet; title ('Co-Occurence matrix');

LCoOc=log2(CoOc);
LCoOc(CoOc==0)=0;  % no nan resulting from inf*0;
H_row=-sum( CoOc.*LCoOc,2 );
H=mean(H_row);

xlabel(strcat('mean entropy for each row is:  ',num2str (H)));
subplot(2,2,4);
modes=sum(CoOc>0,2);
plot(modes);title (strcat('using Thr: ',num2str(Parameter.Spatil.CoOcThr)));
axis([1,size(CoOc,1),0,15]);grid on
xlabel('Labels');ylabel('modes for every cluster')

subplot(2,2,1);
imagesc(Pr_Img);title ('Probabilty to be in Lhat'); colormap jet
xlabel(strcat('mean prob.= ',num2str(mean(Pr)) ));colorbar
subplot(2,2,2)
imagesc(tmp_Labels);title ('Temp Label image'); colormap jet
xlabel(strcat(num2str(length(unique(Lhat))) ,' diffrent labels'))

end