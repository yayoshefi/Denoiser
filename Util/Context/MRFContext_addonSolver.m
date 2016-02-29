function AssignVec2=MRFContext_addonSolver(Data,AssignVec,Centers)
%% AssignVec2=MRFContext_addonSolver(Data,AssignVec,Centers)
% similar to MRFContext but uses a solver from .... complete
%
%
% an mrf solver for minimizing:
%       Data Tem- each Label is given by Dist to Centers
%
%       Smoothness Term- neighboors tend to be similar

global Parameter
K=size(Centers,3);
m=Analysis.LabelsSize(1);
n=Analysis.LabelsSize(2);

[E]=affinity (Data, Centers,0,1);
Dc=exp(reshape(E',[m,n,K]).*10^2);
Sc = ones(K) - eye(K);

[Hc, Vc] = SpatialCues( reshape(Data(1,:),[m,n]) );

gch = GraphCut('open', Dc, Parameter.Spatil.lambda*Sc, exp(-Vc*5), exp(-Hc*5) );
[gch, L] = GraphCut('expand',gch,5);
gch = GraphCut('close', gch);

AssignVec2=L(:)';
end
