function [E_directional,comp]=Compatability (E,AssignVec)
%%  need to update

% Parameters
global Parameter

wsize=Parameter.wsize2^0.5;
row=Parameter.row;      col=Parameter.col;
[K,pnum]=size (E);
m=row-wsize+1;         % number of rows in Assign Mat
n=col-wsize+1;         % number of columns in Assign Mat

if Parameter.spatial.lambda;Parameter.spatial.lambda=0.1; end
% ---- initialization   -------
Lables=col2im(AssignVec,[wsize,wsize],[row,col],'sliding');

% inner pixels (not border)
Logical_C=logical(padarray(ones(m-2,n-2),[1,1]));
E_C=E(ones(K,1),Logical_C(:)');
L_C=reshape( Lables(Logical_C),[m-2,n-2] );

% shifted Right inner matrix
Logical_R=logical( padarray(padarray(ones(m-2,n-2),[1,0]),[0,2],'pre') );
E_R=Parameter.spatial.lambda*E(ones(K,1),Logical_R(:)');
L_R=reshape( Lables(Logical_R),[m-2,n-2] );
comp_R= 2*(L_R==L_C)-1;

% shifted Left inner matrix
Logical_L=logical( padarray(padarray(ones(m-2,n-2),[1,0]),[0,2],'post') );
E_L=Parameter.spatial.lambda*E(ones(K,1),Logical_L(:)');
L_L=reshape( Lables(Logical_L),[m-2,n-2]);
comp_L= 2*(L_L==L_C)-1;

%shifted Up Inner Matrix
Logical_U=logical( padarray(padarray(ones(m-2,n-2),[0,1]),[2,0],'post') );
E_U=Parameter.spatial.lambda*E(ones(K,1),Logical_U(:)');
L_U=reshape( Lables(Logical_U),[m-2,n-2] );
comp_U= 2*(L_U==L_C)-1;

%shifted Down Inner Matrix
Logical_D=logical( padarray(padarray(ones(m-2,n-2),[0,1]),[2,0],'pre') );
E_D=Parameter.spatial.lambda*E(ones(K,1),Logical_D(:)');
L_D=reshape( Lables(Logical_D),[m-2,n-2] );
comp_D= 2*(L_D==L_C)-1;


%{
% E = [C,L,R,U,D] each inner matrix is normalize such sum prob=1
E_Cn=DiagonalMult(E_C,1./sum(E_C,2),'l'); clearvars E_C
E_Ln=DiagonalMult(E_L,1./sum(E_L,2),'l'); clearvars E_L
E_Rn=DiagonalMult(E_R,1./sum(E_R,2),'l'); clearvars E_R
E_Un=DiagonalMult(E_U,1./sum(E_U,2),'l'); clearvars E_U
E_Dn=DiagonalMult(E_D,1./sum(E_D,2),'l'); clearvars E_D
E=[E_Cn,E_Ln,E_Rn,E_Un,E_Dn];
% E=[DiagonalMult(E_C,1./sum(E_C,2),'l'),DiagonalMult(E_L,1./sum(E_L,2),'l'),...
%     DiagonalMult(E_R,1./sum(E_R,2),'l'),DiagonalMult(E_U,1./sum(E_U,2),'l'),...
%     DiagonalMult(E_D,1./sum(E_D,2),'l')];
E(isnan(E))=0;
%}      

E_directional=cat(3,E_R,E_L,E_U,E_D);
comp=cat(3,comp_R(:)',comp_L(:)',comp_U(:)',comp_D(:)');
end