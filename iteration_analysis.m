
h1=figure('name','Iteration Labels');
h2=figure('name','Iteration denoisnig');

K=Parameter.Spatil.MaxIter;
for iter=1:K

    [Context_Output]=removenoise(double(Image),Noise,Analysis.iterations(iter).AssignVec2);
    [r,c]=size(Context_Output);
    Analysis.iterations(iter).result2=psnr(Context_Output,double(Image),255);
    [Analysis.iterations(iter).epsNorm]=...
        ShowCoOc(Analysis.iterations(iter).AssignVec2,false,'EpsNorm','CC');

    wsize=Parameter.wsize2^0.5;
    Lables=col2im(Analysis.iterations(iter).AssignVec2,[wsize,wsize],[r,c],'sliding');
    
    figure(h1);
    subplot(2,round(K/2),iter);
    imshow(Lables,[]);colormap jet;
    title (['iter=',num2str(iter)]);xlabel(['\epsilon-norm',num2str(Analysis.iterations(iter).epsNorm)])
    figure(h2);
    subplot(2,round(K/2),iter);
    imshow(Context_Output,[]);title (['iter=',num2str(iter)]);
    xlabel(['psnr',num2str(Analysis.iterations(iter).result2)])
end
