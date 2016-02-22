function [Centers, Basis,AssignVec,E80]=SplitClusters(Data,Centers,Basis,LargeClusters,AssignVec,E80)
%% [Centers, Basis,AssignVec,E80]=SplitClusters(Centers,Basis,LargeClusters,MedianProjection,AssignVec,E80,normalize)
% splits each cluster center on MainVec median
% the half smaller than MEDIAN replaces index
% the half larger than MEDIAN is an addtive index

global Parameter
[wsize2,Dim,K]=size(Basis);
if Parameter.normalize==2
    Mvec=2;
else Mvec=1;
end
cost=5000;

LargeClusters=sort(LargeClusters,'descend');
Dbugger=zeros(size(LargeClusters));

while ~isempty(LargeClusters)

    Curr=LargeClusters(1);
%%%
    clstindx=find((AssignVec==Curr));
    clstsize=sum(AssignVec==Curr);
    cluster_0=Data(:,clstindx)-Centers(:,ones(1,clstsize),Curr);
    
    switch Parameter.SplitType
        case 'median'
            projectionOnMainVec=Basis(:,Mvec,Curr)'*cluster_0;
            clusterMed=median(projectionOnMainVec); 
            
            clst2=projectionOnMainVec>clusterMed;
            clst1=~clst2;
            size2=sum(clst2);       size1=clstsize-size2;
            if size2<Parameter.minclstsize/2
                LargeClusters(1)=[];continue
            end
            
            NewCenter2=mean(Data(:,clstindx(clst2)),2);
            NewCenter1=mean(Data(:,clstindx(clst1)),2);
            
        case 'totvar'
%             opts = statset('Display','final','MaxIter',100);
%             [indx, C]=kmeans(Data(:,clstindx)',2,'Options',opts);           
            [indx, C]=kmeans(Data(:,clstindx)',2);
            clst1=indx==1;          clst2=indx==2;
            size1=sum(indx==1);     size2=sum(indx==2);
            if size1<Parameter.minclstsize || size2<Parameter.minclstsize  %|| sum(TotVar)<cost
                LargeClusters(1)=[];continue;
            end
            NewCenter2=C(1,:)';
            NewCenter1=C(2,:)';
            
    end

    clst1_0= Data(:,clstindx(clst1))-NewCenter1(:,ones(1,size1));
    clst2_0= Data(:,clstindx(clst2))-NewCenter2(:,ones(1,size2));

    [U1,S1,~]=svd(clst1_0,'econ'); 
    [U2,S2,~]=svd(clst2_0,'econ');
    
    NewBasis1=zeros(wsize2); NewBasis2=zeros(wsize2);
    NewBasis1(:,1:size(U1,2))=U1;
    NewBasis2(:,1:size(U2,2))=U2;
            
    if strcmp(Parameter.SplitType,'totvar')
        TotVar=mean( (Basis(:,1:wsize2,Curr)'*cluster_0).^2 ,2 );
%         withinVar=( ( mean((NewBasis1'*clst1_0).^2,2) )*size1+( mean((NewBasis2'*clst2_0).^2,2) )*size2 )/clstsize;
        withinVar= sum( mean((U1'*clst1_0).^2,2)*size1/clstsize)+...
            sum( mean((U2'*clst2_0).^2,2)*size2/clstsize);
        betweenVar=sum(TotVar)-sum(withinVar);
        
        Dbugger(numel(LargeClusters))=betweenVar;
        

        if betweenVar<cost
            LargeClusters(1)=[];
            continue;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Var(P)=E[Var(P|Centers)]+Var(E[P|Centers])  == within + between
%between= TotalVar-within
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    AssignVec(clstindx(clst2))=K+1;
    K=K+1;

    Centers=cat(3,Centers,NewCenter2);
    Centers(:,:,Curr)=NewCenter1;

    Basis=cat (3,Basis,NewBasis2);
    Basis(:,:,Curr)=NewBasis1;                

    S2=diag(S2).^2;         E2=cumsum(S2);
    S1=diag(S1).^2;         E1=cumsum(S1);
    Eper1=E1./E1(end);      Eper2=E2./E2(end);
    E80_1=find(Eper1>0.8,1);
    E80_2=find(Eper2>0.8,1);
    E80=cat(3,E80,E80_2);
    E80(Curr)=E80_1;

    LargeClusters(1)=[];

end
if ~isempty(Dbugger)&&strcmp(Parameter.SplitType,'totvar')
    disp(['betweenVar threshold is: ',num2str(cost)])
    disp (num2str(round(Dbugger)))
end

end
