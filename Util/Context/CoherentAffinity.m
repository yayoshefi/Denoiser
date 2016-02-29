function [E]=CoherentAffinity (Data, Centers,varargin)
%% [E]=affinity (Data, Centers,ImageSize,varargin)
% calculates the affinity matrix for coherent clustring Y. Hel-Or 6/2014
% Data is [wsize2 X n] the whole image patches.
% Centers is [wzise2 X K] each cluster mean.
%
% The affinity matrix E is a vercat of 5 matricies, each one is  normalized
% such the sum of prob. for a random walker to get to a Cluster equal 1.
%
% parmeters for the distance computing is enterd in Pairs :
% {'metric'  'dim' 's' 'basis' 'lambda'  'sigma'  'normalize'}
% currently works for only the 'euclidean' case. dim=0
% lambda sets the ratio between Binary and Unitary
% FUTURE UPDATES:
% set lambda prameter for choosing Binary Vs Unitary dist ratio
% set sigma for the energy dist function
% enter rows and columns in the input


% Parameters
global Parameter
[wsize2, pnum]=size(Data);
K=size(Centers,3);

row=Parameter.row;      col=Parameter.col;

m=row-wsize2^0.5+1;         % number of rows in Assign Mat
n=col-wsize2^0.5+1;         % number of columns in Assign Mat

pnames = {  'dim'        's'         'basis' 'lambda',      'sigma' };
dflts =  { zeros(1,K), ones(wsize2,1,K),0,   0.3,        wsize2^(1/4)};
[metric,Dim,S,Basis,lambda,sigma,normalize] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});
if Parameter.normalize==0
    sigma=3.5*sqrt(wsize2);end
    %similarity parameter for the gaussian similarity
% sigma=wsize2/1.41;

% initialization
% % E_C=zeros((m-2)*(n-2),K);E_R=zeros((m-2)*(n-2),K);E_L=zeros((m-2)*(n-2),K);
% % E_U=zeros((m-2)*(n-2),K);E_D=zeros((m-2)*(n-2),K);
% % Thr=0.3;
knn=0; value=true;
[E]=affinity (Data, Centers,knn,value);          %,varargin{1:6});
LogicalInner_C=logical(padarray(ones(m-2,n-2),[1,1]));
E_C=E(LogicalInner_C);
%shifted Right inner matrix
LogicalInner_R=logical( padarray(padarray(ones(m-2,n-2),[1,0]),[0,2],'pre') );
E_L=lambda*E(LogicalInner_R);
%shifted Left inner matrix
LogicalInner_L=logical( padarray(padarray(ones(m-2,n-2),[1,0]),[0,2],'post') );
E_R=lambda*E(LogicalInner_L);
%shifted Down Inner Matrix
LogicalInner_D=logical( padarray(padarray(ones(m-2,n-2),[0,1]),[2,0],'pre') );
E_U=lambda*E(LogicalInner_D);
%shifted Up Inner Matrix
LogicalInner_U=logical( padarray(padarray(ones(m-2,n-2),[0,1]),[2,0],'post') );
E_D=lambda*E(LogicalInner_U);

% E = [C,L,R,U,D] each inner matrix is normalize such sum prob=1
E_Cn=DiagonalMult(E_C,1./sum(E_C,2),'l'); clearvars E_C
E_Ln=DiagonalMult(E_L,1./sum(E_L,2),'l'); clearvars E_L
E_Rn=DiagonalMult(E_R,1./sum(E_R,2),'l'); clearvars E_R
E_Un=DiagonalMult(E_U,1./sum(E_U,2),'l'); clearvars E_U
E_Dn=DiagonalMult(E_D,1./sum(E_D,2),'l'); clearvars E_D
E=[E_Cn,E_Ln,E_Rn,E_Un,E_Dn];
% E=[DiagonalMult(E_C,1./sum(E_C,2),'l'),DiagonalMult(E_L,1./sum(E_L,2),'l'),...
%     DiagonalMult(E_R,1./sum(E_R,2),'l'),DiagonalMult(E_U,1./sum(E_U,2),'l'),...
%     DiagonalMult(E_D,1./sum(E_D,2),'l')];
E(isnan(E))=0;
        

% E=sparse(E);
end