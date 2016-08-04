function [AssignVec2]=MRFContext(Data,AssignVec,Centers)

global Parameter
row=Parameter.row;      col=Parameter.col;
wsize=sqrt(Parameter.wsize2);
K=max(AssignVec);

AssignImg=col2im(AssignVec,[wsize,wsize],[row,col],'sliding');
[m,n]=size(AssignImg);

%calc evidence
[E]=affinity (Data, Centers,0,1);

%calc compatibilty: co-occurence
glcm = graycomatrix(AssignImg,'GrayLimits',[1,K],'NumLevels',K,'Offset',[0,1;1,0]);
% glcm (i,j,1)= #lablel j to the right of label i (each row)
% glcm (i,j,1)= #lablel i to the left of label j (each column)
% the glcm (i,j,2) is: rows, # lables below (+1 to y spatial subscript) 

Prob_RN=normalizecomatrix(glcm(:,:,1),'row');
Prob_DN=normalizecomatrix(glcm(:,:,2),'row');
Prob_LN=normalizecomatrix(glcm(:,:,1),'col');
Prob_UN=normalizecomatrix(glcm(:,:,2),'col'); % probabily for your up neigbour


% [s,t,weights,p]=GraphWeights (AssignVec); %early obsolete version


    
AssignVec2=AssignVec;
p=nodesIndx(m,n); %Sourse node=1, Target(sink) node= m*n
SourceUnitary=E(1,:);
tic;
h=waitbar(0,'bulinding Graph & Finding min cut..','CreateCancelBtn','setappdata(gcbf,''cancel'',1)');
setappdata(h,'cancel',0);

% aplha-expension graph cut one loop
for label=2:K
Source=1;       Target=label;   %Labels for S,T

    if getappdata(h,'cancel')
    break
    end
    waitbar(label/K);
    
s=[ones(1,(m-2)*(n-2))  ,repmat(p,1,4)          ,p];
% Source X pixels       pixels X 4 conectivty   pixels
t=[p                    ,[p+n,p+1,p-n,p-1]      ,(m*n)*ones(1,(m-2)*(n-2))];
% pixels                neigbours               Target X pixels
% evidence to Source    compatabilty            evidence to Target

right=sub2ind([K,K],AssignVec(p),AssignVec(p+n));
down=sub2ind([K,K],AssignVec(p),AssignVec(p+1));
ToTheLeft=sub2ind([K,K],AssignVec(p),AssignVec(p-n));
up=sub2ind([K,K],AssignVec(p),AssignVec(p-1));

compR=Parameter.spatial.lambda*Prob_RN(right); compD=Parameter.spatial.lambda*Prob_DN(down);
compL=Parameter.spatial.lambda*Prob_LN(ToTheLeft); compU=Parameter.spatial.lambda*Prob_UN(up);

weights=[SourceUnitary(p)    ,compR,compD,compL,compU    ,E(Target,p)];
G=digraph(s,t,weights); %undefind for matlab 2014


[mf,~,cs,ct] = maxflow(G,1,m*n);


changedNodes=ct(  ct>m+1 & ct <=( m*(n-1) )  );  %in the ct half & within nodes
changedNodes(mod (changedNodes,m)==0)=[]; 
changedNodes(mod(changedNodes-1,m)==0)=[]; %delete 1st and last row nodes

AssignVec2 (changedNodes)=Target;
% new=sub2ind(size(E),Target(ones(size(changedNodes))),changedNodes);
SourceUnitary(changedNodes)=E(2,changedNodes);
%  the new Assignment is for all points change thier current assignment

end
delete(h);
disp(['one loop time: ',num2str(toc)])
% for evry iteration if changes->update AssignVec2
%build the new comp functions
end

function p=nodesIndx(m,n)
% nodes index list, for pixel location in mXn inner matrix
p=1:m*(n-1);
p(1:m)=[];
p(mod (p,m)==0)=[];
p(mod(p-1,m)==0)=[];

end

function glcmN=normalizecomatrix(glcm,direction)

switch direction
    case 'row'
        glcmnorm=sum(glcm,2);
        glcmN=glcm./glcmnorm(:,ones(1,size(glcm,2)));
    case 'col'
        glcmnorm=sum(glcm,1);
        glcmN=glcm./glcmnorm(ones(size(glcm,1),1),:);
end       
end
