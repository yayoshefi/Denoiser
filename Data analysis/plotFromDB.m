fprintf('Plotting statistics for %i clusters and noise std=%i\n',full_Data.Parameters.values.kmeans,full_Data.Parameters.sigma)
% PSNR [K-means, CoC , ORACLE]'
Psnr=struct2cell(full_Data.Psnr);
Psnr= squeeze( cell2mat(Psnr(2:4,:,:)) );

% Sparsity [K-means, CoC , ORACLE]'
Sparsity=struct2cell(full_Data.Sparsity);
Sparsity= squeeze( cell2mat(Sparsity([3,7,5],:,:)) );


%% sparsity Vs Psnr per each image
K=full_Data.Parameters.values.kmeans^2;
figure;title('Psnr Vs Sparsity')
plot( Sparsity(:,1:end-1)/K,Psnr(:,1:end-1) );
xlabel('Sparse'); ylabel ('Psnr')

%% CoC change, mu,sig 
improve=Psnr(2,:)-Psnr(1,:);

sig=std(improve(1:end-1)); mn=mean(improve(1:end-1));
figure; title('change between K-means and CoC')
stem (improve(1:end-1));
xlabel({'image number',['\mu=',num2str(mn,'%1.4f'),' ; \sigma=',num2str(sig,'%1.4f')]})

%% bar plot of mean sparsity and psnr (instead of latex)
x={'K-means','CoC','ORACLE'};

figure;
[AX,H1,H2] =plotyy( 1:3:9,Sparsity(:,end),2:3:9,Psnr(:,end),'bar','bar' );
H1.BarWidth=0.3;                            H2.BarWidth=0.3;
AX(1).YLim;                                 AX(2).YLim=[22,29];
AX(1).YTick;                                AX(2).YTick=22:2:30;
c1=AX(1).ColorOrder(1,:);                   c2=AX(1).ColorOrder(2,:);
H1.FaceColor=c1;                            H2.FaceColor=c2;
ylabel(AX(1),['|| ',full_Data.Parameters.CoOc.Type,' ||_{\epsilon}'],'Color',c1);    ylabel(AX(2),'Psnr','Color',c2);
legend({'Sparsity','Psnr'})

AX(1).XTickLabel=x;

 text(1:3:9,Sparsity(:,end),num2str(Sparsity(:,end)/K,'%0.2f'),'HorizontalAlignment','center',... 
'VerticalAlignment','bottom')

%% worst and best improvments
[value, index]=sort(improve);
disp ('worst images')
disp ([full_Data.Psnr(index(1)).Name,'  ',full_Data.Psnr(index(2)).Name,'  ',full_Data.Psnr(index(3)).Name])

disp ('best images')
disp ([full_Data.Psnr(index(end)).Name,'  ',full_Data.Psnr(index(end-1)).Name,'  ',full_Data.Psnr(index(end-2)).Name])

