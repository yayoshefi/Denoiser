function AssignVec2=SpectralContext (Data,AssignVec,Centers,varargin)
%% AssignVec2=SpectralContext (Data,AssignVec,Centers,varargin)
%  the function uses Spectral clustering with Neighbors affinity matrix
%  such that there are 5 affinitys (C,R,L,U,D) from pixels to Labels and
%  the Spectral clustering can use all of them in order to find the
%  shortest distance in the emmbeded space

global Parameter

Parameter.metric ='euclidean';
K=size(Centers,3);
wsize=sqrt(Parameter.wsize2);
m=Analysis.LabelsSize(1);
n=Analysis.LabelsSize(2);


pnames = {'s'           'basis' };
dflts =  {ones(wsize^2,1,K),0 };
[Energy,Basis] = internal.stats.parseArgs(pnames, dflts, varargin{:});

%if Parameter.Spatil.lambda==0; Parameter.Spatil.lambda=0.1; end
%E=CoherentAffinity (Patches, Centers,'s',Energy,'basis',Basis,'lambda',lambda);
[E]=affinity (Data, Centers,20,1);
B = Spatialaffinity ();
E_spatial=[E;Parameter.Spatil.lambda*B];
[AssignVec2]=MySpectralClustering(E_spatial,'e');
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
