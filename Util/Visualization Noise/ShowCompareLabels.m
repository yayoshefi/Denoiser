function [Labeldist]=ShowCompareLabels(L,Data,Oracle)
global Parameter Analysis
if isfield(Analysis,'ORACLE');       Oracle=Analysis.ORACLE;end
   
row=Parameter.row;      col=Parameter.col;      wsize=Parameter.wsize2^0.5;
if isvector(Oracle)
    Oracle = col2im(Oracle,[wsize wsize],[row col],'sliding');end
if isvector(L)
    L = col2im(L,[wsize wsize],[row col],'sliding');end

ne=(Oracle~=L);
K=max(Oracle(:));

[Center1,~,Lhat1,~,~]=UpdateCenter(Data,Oracle(:),false);
[Center2,~,Lhat2,~,~]=UpdateCenter(Data,L(:),false);
[order, D]=Dist2SubSpace (squeeze(Center2),Center1,'knn',K);
% D(M,J)- is the dist between Center M in L and Center J in the Oracle

ind=sub2ind([K,K],Oracle(:),L(:));
Labeldiff=D(ind);

Labeldiff=col2im(Labeldiff,[wsize wsize],[row col],'sliding');
% the distance is based on the ground distance between center is the noisy
% image, and not from the clean image centers.

Labeldist = mean(Labeldiff(:));
figure;
subplot (1,2,1)
imshow(ne);title ('missplaced labels')
subplot (1,2,2)
imshow(log(Labeldiff),[]);title ('distance between labels')
xlabel(['computed distance: ',num2str(Labeldist)])
end