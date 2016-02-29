function AssignVec2= MinEntropyCC (Data,AssignVec,Centers)
%%AssigVec2= MinEntropyCC (Data,AssignVec,Centers)
% the function minimize the Co-Occurence sparsity  by trying to converge 
% to delta function to as an optimum.
% this works only at hist space without changing Visual space

global Parameter

wsize=sqrt(Parameter.wsize2);
K=size(Centers,3);
AssignImg=col2im(AssignVec,[wsize,wsize],[Parameter.row,Parameter.col]);
NN=Parameter.Spatil.NN;
Neigbour=im2col(AssignImg,[NN,NN],'sliding');
Neigbour(ceil(NN^2/2),:)=[];
LocalHist=histc(Neigbour,1:K,1);        

Histoids=K*eye(K);
[AssignVec2,C]=kmeans(LocalHist',K,'start',Histoids,'distance','correlation');
end