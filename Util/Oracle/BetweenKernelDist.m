function [Dist2R,Dist2D,Dist2L,Dist2U]=BetweenKernelDist (AssignVec) 
%% [Dist2R,Dist2D,Dist2L,Dist2U]=BetweenKernelDist (AssignVec)
% This is the Compatable function between each cluster with direction.
% AssignMat can be either a matrix or even a vector
% IF AssignMat is a vector, entering dim and wsize is required
% dim = [row,col], or the output of size (Image)
% wsize is the window size used to open the window to patches

global Parameter

row=Parameter.row;      col=Parameter.col;
wsize=sqrt(Parameter.wsize2);
AssignMat=col2im(AssignVec,[wsize,wsize],[row,col],'sliding');

K=max(AssignVec(:));
[m, n]=size(AssignMat);

Dist2R=zeros(K);Dist2D=zeros(K);Dist2L=zeros(K);Dist2U=zeros(K);

%too svere, this delte all border pixels
p=1:m*(n-1);                % index of all Image except last column
p(1:m)=[];                  % delete first column indexes   
p(mod (p,m)==0)=[];         % delete bottom row
p(mod(p-1,m)==0)=[];        % delele top row

clst2R=AssignMat(p+m);      %what cluster does your right neighbour belong to
clst2D=AssignMat(p+1);
clst2L=AssignMat(p-m);
clst2U=AssignMat(p-1);

pixel_num=(m-1)*(n-1);

for k=1:K
    Curr=(AssignMat(p)==k);
    Dist2R(k,:)=histc(clst2R(Curr),1:K)/sum(Curr);
    Dist2D(k,:)=histc(clst2D(Curr),1:K)/sum(Curr);
    Dist2L(k,:)=histc(clst2L(Curr),1:K)/sum(Curr);
    Dist2U(k,:)=histc(clst2U(Curr),1:K)/sum(Curr);  %this is the ONLY conditianal prob. the rest is normalizes to the whole amount
    %/pixel_num - maybe insted of /sum(Curr)
end

end