function [psnrMat]=PlotPsnr_Center(Method,REFMethod,wsizes)
%% []=PlotData(psnr,RefPsnr,wsizes)
% takes struct psnr and plots it to graph

W=length(wsizes);
psnrMat=[];
figure; 
for w=1:W
    
    [Methodpsnr,MethodcentersCount]=struct2data (Method);
    
    [REFpsnr,REFcentersCount]=struct2data (REFMethod);
    low=min([MethodcentersCount,REFcentersCount])-5;
    up=max([MethodcentersCount,REFcentersCount])+5;

    subplot(W,1,w); hold on;
    
    [Hax, H1,H2]=plotyy(10:10:10*length(Methodpsnr),Methodpsnr,10:10:10*length(Methodpsnr),MethodcentersCount);
    set(H1,'LineStyle','-','Color','b','LineWidth',2);
    set (H2,'LineStyle',':','Color','b','LineWidth',0.2) 
    grid (Hax(1),'on');
    set(Hax(1),'Color','w')
    set(Hax(2),'YColor','k','YLim',[low,up],'YTick',linspace(low,up,5));
    
    [Hax, H1,H2]=plotyy(10*[1:length(REFpsnr)],REFpsnr,10*[1:length(REFpsnr)],REFcentersCount');
    set(H1,'LineStyle','-.','Color','r','LineWidth',2);
    set (H2,'Color','r','LineStyle',':','LineWidth',0.2);
    set(Hax(2),'YColor','k','YLim',[low,up],'YTick',linspace(low,up,5));

    
    ylabel(Hax(1),['window size=',num2str(wsizes(w))])
    

    psnrMat=cat(1,psnrMat,Methodpsnr);

end
xlabel('sigma');
subplot(W,1,1);legend('Method','Reference')
ylabel('psnr')


function [psnr,centersCount]=struct2data (name)
cell_psnr=struct2cell( name.psnr.(['wsize',num2str(wsizes(w))]) );
psnr=cell2mat(cell_psnr)';
cell_centers=struct2cell( name.centers.(['wsize',num2str(wsizes(w))]) );
centersCount=cell2mat(cell_centers)';

end

end

