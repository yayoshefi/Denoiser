function AssignVec2=MutualDist(Data,AssignVec,Centers)
%% AssignVec2=MutualDist(Data,AssignVec,Centers)
% this function returns assignment that is a result of minimizing the
% mutual distance both in visul space and both in histogrma space.
% 
% in visual space we use a simple euclidean dist (L2 norm)

% in Hist spcae we use Earth Mover's Dist based on L2 norm between the
% features vectors (represented as Centers). the signiture is the histogram
% of the neighborhood based on distance between Centers


global Parameter Analysis

K=size(Centers,3);
wsize=sqrt(Parameter.wsize2);
m=Analysis.LabelsSize(1);
n=Analysis.LabelsSize(2);


[dim,pnum]=size (Data);
NN=Parameter.spatial.NN;  %window2
padding=floor(NN/2);

Lhat=AssignVec; MaxIter=3; figure;
for iter=1:MaxIter
    AssignImg=col2im(Lhat,[wsize,wsize],[Parameter.row,Parameter.col]);
    
    if Analysis.DebuggerMode && ~(mod(iter+1,1));  % Debug
        subplot(MaxIter,1,iter);imagesc(AssignImg);colormap (Analysis.ColorMap);
        axis image off ;title(['iteration ',num2str(iter)]);
    end
    AssignImg=padarray(AssignImg,[padding,padding],-1);
    
    Neigbour=im2col(AssignImg,[NN,NN],'sliding');
    Neigbour(ceil(NN^2/2),:)=[];
    Hist=histc(Neigbour,1:K,1)'/(NN^2-1);
%     p=logical(padarray(ones(m-2*padding,n-2*padding),[padding,padding]));
%     [C,H]=clusterRep(Data,Lhat(p),Hist);
    [C,H]=clusterRep(Data,Lhat,Hist);

    MutualDist=inf*ones(1,pnum);  %AssignVec2=Lhat;
    for k=1:K
        d_vis=sqrt(  abs( sum(Data.^2)-2*C(:,:,k)'*Data+C(:,:,k)'*C(:,:,k) )  );
        d_hist=FastEMD(squeeze(C),H(k,:),squeeze(C),Hist);
        
        tempDist=d_vis+Parameter.spatial.lambda*d_hist;
        index=tempDist<MutualDist;
        MutualDist(index)=tempDist(index);
        Lhat(index)=k;
    end
end
AssignVec2=Lhat;
end
