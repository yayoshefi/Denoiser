function [Labeldist]=ShowCompareLabels(L1,L2,Data)
global Parameter Analysis
if isfield(Analysis,'ORACLE');       L2=Analysis.ORACLE;end
   
row=Parameter.row;      col=Parameter.col;      wsize=Parameter.wsize2^0.5;
if isvector(L2)
    L2 = col2im(L2,[wsize wsize],[row col],'sliding');end
if isvector(L1)
    L1 = col2im(L1,[wsize wsize],[row col],'sliding');end

ne=(L2~=L1);
K=max([L1(:);L2(:)]);

if exist('Data','var')
    [Center1,~,Lhat1,~,~]=UpdateCenter(Data,L1(:),false);
    [Center2,~,Lhat2,~,~]=UpdateCenter(Data,L2(:),false);
    [order, D]=Dist2SubSpace (squeeze(Center2),Center1,'knn',K);
    % D(M,J)- is the dist between Center M in L and Center J in the Oracle

    ind=sub2ind([K,K],L2(:),L1(:));
    Labeldiff=D(ind);

    Labeldiff=col2im(Labeldiff,[wsize wsize],[row col],'sliding');
    % the distance is based on the ground distance between center is the noisy
    % image, and not from the clean image centers.

    Labeldist = mean(Labeldiff(:));
end

figure;
subplot (1,2,1)
imshow(ne);title ('missplaced labels')

if exist('Data','var')
    subplot (1,2,2)
    imshow(log(Labeldiff),[]);title ('distance between labels')
    xlabel(['computed distance: ',num2str(Labeldist)])
end
end