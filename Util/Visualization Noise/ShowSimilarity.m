function []=ShowSimilarity (E,AssignVec)
%% []=ShowSimilarity (E)


[K, pnum]=size (E);
if length(AssignVec)~=pnum
    error ('myApp:argChk','Assign Vector and affinity are not same length')
end
perm=randperm (pnum,round(sqrt(pnum)));
e=E(:,perm);
[~,I]=sort (AssignVec(perm));
Sim= e(:,I)'*e(:,I);
figure; imagesc(Sim)
end