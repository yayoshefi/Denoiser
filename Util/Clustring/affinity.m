function [E]=affinity (Data, Centers,knn,value,varargin)
%% [E]=affinity (Data, Centers,knn,value,varargin)
%
% the affinity matrix E is [K by pnum], inconsistent with Toky (directional)
% calculates the affinity matrix. the affinity is  normalized to prob.1
% Data is [wsize2 X n] the whole image patches.
% Centers is [wzise2 X 1 X K] each cluster mean.
%
% Knn sets how much non empty elements in each col.
% if knn=0 then a complete Affinity is computed between all Centers.
%
% value sets if the affinity gets the gaussian func value or just a logical
%
% parmeters for the distance computing is enterd in Pairs :{'metric'  'dim'
% 's' 'basis'}
% currently works for only the 'euclidean' case. dim=0
% FUTURE UPDATES:
% set sigma for the energy dist function
% enter rows and columns in the input


% Parameters
global Parameter
[wsize2,pnum]=size (Data);
[K]=size(Centers,3);
row=Parameter.row;      col=Parameter.col;

m=row-wsize2^0.5+1;         % number of rows in Assign Mat
n=col-wsize2^0.5+1;         % number of columns in Assign Mat

% E=zeros((m-2)*(n-2),K);
E=zeros(K,pnum);
if Parameter.normalize==2
    sigma=wsize2/4;
else    sigma=wsize2/4;
end

pnames = {'dim'       's'   'basis'};
dflts =  {zeros(1,K), ones(wsize2,1,K),0};
[Dim,S,Basis] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});


if knn==0; knn=K;
else knn=min(knn,K); end
S=PseudoDiagInv(S);

[AssignVec, Distances]=Dist2SubSpace (Data,Centers,'knn',knn,varargin{:});     %case basis is in varargin can be a problem
col_sub=meshgrid(1:pnum,1:knn);
ind=sub2ind([K,pnum],AssignVec(:),col_sub(:));
if value==0
    E(ind)=1;
else
    if strcmp(Parameter.metric,'euclidean')
        E(ind)=exp(Distances(:)/(-sigma^2));
    else E(ind)=Distances(:);           %mahalanobis
    end
end
cumE=sum(E,1);
NormalizeFactor=1./cumE';
NormalizeFactor( cumE'==0 )=0;                  %outlier point- far from all centers
if sum (cumE==0)>log(pnum)
    error ('AffinityFunc:Outliers','%g points where Sum[Pr(k)]=0 out of %g total points',sum (cumE==0),pnum)
end
E=DiagonalMult(E,1.*NormalizeFactor,'r');           %normalize the affinity matrix


end