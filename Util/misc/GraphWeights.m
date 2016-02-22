function [S,T,Weights,p]=GraphWeights (AssignVec)
%% [S,T,Weights,p]=GraphWeights (AssignVec)
% computes the Graph components
% using Dist2R....


global Parameter

m=Parameter.row-sqrt(Parameter.wsize2)+1; n=Parameter.col-sqrt(Parameter.wsize2)+1;
K=max(AssignVec);

[Dist2R,Dist2D,Dist2L,Dist2U]=BetweenKernelDist (AssignVec) ;

% Source index list
p=1:m*(n-1);
p(1:m)=[];
p(mod (p,m)==0)=[];
p(mod(p-1,m)==0)=[];

% Source and Target for each direction
S_R=p;
T_R=p+n;
index_R=sub2ind([K,K],AssignVec(S_R),AssignVec(T_R));
W_R=Dist2R(index_R);

S_D=p;
T_D=p+1;
index_D=sub2ind([K,K],AssignVec(S_D),AssignVec(T_D));
W_D=Dist2D(index_D);

S_L=p;
T_L=p-n;
index_L=sub2ind([K,K],AssignVec(S_L),AssignVec(T_L));
W_L=Dist2L(index_L);

S_U=p;
T_U=p-1;
index_U=sub2ind([K,K],AssignVec(S_U),AssignVec(T_U));
W_U=Dist2U(index_U);

% combine all
Weights=[W_R;W_D;W_L;W_U];
Weights=Weights(:)';

T=[T_R;T_D;T_L;T_U];
T=T(:)';

S=[p;p;p;p];
S=S(:)';
end