%% checking nystrom method
% check if the nystrom methid is a good approximation to regular k means
% use after Patches initalization stage clustering


perm=randperm (length(Patches));
[~,reorder]=sort(perm);
samp=round(length(Patches)/150);
% samp=1e3;
LandMarks=Patches(:,perm(1:samp));
E=zeros(samp,length(Patches)-samp);
band=500;
for bp=1:band:samp
    chunk=bp:min(bp+band-1,samp);
    E(chunk,:)=affinity (Patches(:,samp+1:end),reshape(LandMarks(:,chunk),[],1,length(chunk)),0,1);   
end

[AssignVecNys]=MySpectralClustering(E,'nystrom','Centers',LandMarks);
AssignVecNys=AssignVecNys(reorder);
[Nystrom_Output]=removenoise(double(Image),Noise,AssignVecNys);

resultNys=psnr(Nystrom_Output,double(Image),255);
Out={Output,Nystrom_Output};
res={result,resultNys};

ShowClusters(Image,Noise,AssignVecNys,Centers,Energy);  %Centers & Energy are not updated
PrintDnoise (Out,res);