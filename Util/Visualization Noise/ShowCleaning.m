function ShowCleaning (Image,Noise,Clean)
%% [handle]=ShowCleaning (Image,Noise,Clean)
%
% handle is the figure handle for the shown figure
% Image is the input image, Noise is the currently random noise
% Clean is the patch estimate after shrinkage (and before accumaliting) 

global Parameter Analysis
wsize=Parameter.wsize2^0.5;

bins=1000;

[VarP,PatchSNR]=ShowPatchSNR (Image,Noise);
Analysis.EmpiricalPathSNR=PatchSNR;
Analysis.CleanSignalVariance=VarP;

P=im2col(Image,[wsize,wsize],'sliding');
N=im2col(Noise,[wsize,wsize],'sliding');

PatchRmse= sqrt( mean((P-Clean).^2) );

Analysis.DC_Change=mean(P)+mean(N)-mean(Clean);     %DC of noisy signal - DC after Clean

% -------- ploting --------------
UVar=max(VarP); LVar=min(VarP)-0.01;
VarScale=(UVar-LVar)/(bins+1);

USNR=max(PatchSNR); LSNR=min(PatchSNR)-0.01;
SNRscale=(USNR-LSNR)/(bins+1);

X_Var=zeros(1,bins);            X_SNR=zeros(1,bins);
Rmse_Var=zeros(1,bins);         Rmse_SNR=zeros(1,bins);
for i=1:(bins+1)
    X_Var(i)=mean(VarP(  (LVar+(i-1)*VarScale)< VarP & VarP <=(LVar+(i)*VarScale)  ));
    Rmse_Var(i)=mean(PatchRmse(  (LVar+(i-1)*VarScale)< VarP & VarP <=(LVar+(i)*VarScale)  ));
    
    X_SNR(i)=mean(PatchSNR(  (LSNR+(i-1)*SNRscale)< PatchSNR & PatchSNR <=(LSNR+(i)*SNRscale)  ));
    Rmse_SNR(i)=mean(PatchRmse(  (LSNR+(i-1)*SNRscale)< PatchSNR & PatchSNR <=(LSNR+(i)*SNRscale)  ));    
end
C=cov([PatchRmse(:),VarP(:),PatchSNR(:)]);
disp (['cov[Rmse,signal Var]=',num2str(C(1,2))]);
disp (['cov[Rmse,Patch SNR]=',num2str(C(1,3))]);

Analysis.Handles.Cleaning=figure;

subplot (2,2,1)
plot (X_Var,Rmse_Var)
xlabel('Var(Patch)');           ylabel ('Rmse')

subplot (2,2,3)
plot (X_SNR,Rmse_SNR)
line ([0.45,0.45],[0,max(Rmse_SNR)],'Color',[1,0,0]);
xlabel('Patch SNR');           ylabel ('Rmse')

end