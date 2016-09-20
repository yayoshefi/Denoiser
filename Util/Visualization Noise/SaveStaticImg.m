function []=SaveStaticImg(ORACLE,ORACLEOutput,ORACLEresult,BM3dresult, BM3dOutput,...
                            name,sigma,InputsImg,InputsStr)
global Parameter
wsize=sqrt(Parameter.wsize2);
row =Parameter.row; col=Parameter.col;
src=strcat(Parameter.location,'\Results\',date,'\',name,'_sigma',num2str(sigma),'\static');
mkdir(src);

figure;imshow(col2im(ORACLE,[wsize,wsize],[row,col]),[]);colormap jet;title (['ORACLE for w_1=',num2str(wsize)])
xlabel(['Psnr ',num2str(ORACLEresult)]);
saveas(gcf,strcat(src,'\ORACLE.png'))

CoOc_V1 (lcm(ORACLE,Parameter.spatial.AssginType,Parameter.spatial.CoOc),true);
saveas (gcf,strcat(src,'\ORACLECoOc.png'))

figure; imshow(ORACLEOutput,[]); title ('ORACLE Denoising')
saveas(gcf,strcat(src,'\ORACLEdenosing.png'))

if exist('BM3dOutput','var')
figure; imshow(BM3dOutput,[]); title (['BM3D Denoising (',num2str(BM3dresult),')'])
saveas(gcf,strcat(src,'\DM3Ddenosing.png'))
end

PrintDnoise (InputsImg,InputsStr)


trgt=strcat('C:\Users\Yair\Dropbox\Thesis\Latex\figures\',date,'\',name,'_',Parameter.spatial.CoOc,'_sigma',num2str(sigma));
copyfile(src,trgt)

end