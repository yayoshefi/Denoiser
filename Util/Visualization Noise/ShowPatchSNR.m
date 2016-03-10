function [VarP,PatchSNR]=ShowPatchSNR (Image,Noise)
%% [VarP,PatchSNR]=ShowPatchSNR (Image,Noise)
% calculate the empirical variance
% calculate the empirical SNR for each patch


global Parameter Analysis
wsize=sqrt(Parameter.wsize2);
Analysis.P=im2col(double(Image),[wsize,wsize],'sliding');
Analysis.N=im2col(double(Noise),[wsize,wsize],'sliding');

[~,pnum]=size(Analysis.P);
if Parameter.wsize2==1;     VarP=Analysis.P.^2;
else                        VarP=var(Analysis.P,1);
end
if Parameter.wsize2==1;     VarN=Analysis.N.^2;
else                        VarN=var(Analysis.N,1);
end


PatchSNR=sqrt(VarP./VarN);

end