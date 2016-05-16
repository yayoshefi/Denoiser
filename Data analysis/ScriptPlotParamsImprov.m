%% script for making Graphs to learn parameters
global Analysis
clearvars data_lambda data_NN data_reference
% change if necessary: struct name,strctfields,strctvalues

%%
Params=Analysis.Arrays;
% Params=Arrays;
%                   use Arrays- when only 1 parameter need to compare, edit
%                   the Array sructure and not the original
ToPlot=1:length(Params.sigma);
%%
Params.Lambda=1:length(Params.Lambda);
Params.CoOcThr=1:length(Params.CoOcThr);
sigma_array=Params.sigma;
Params=rmfield(Params,'sigma');
fields = fieldnames(Params);

%non relevent
Thr_index=1;
Params.CoOcThr=Params.CoOcThr(Thr_index);   Params.normalize=Params.normalize(1);


for w=1:length(Params.wsize)
    Params.wsize=Arrays.wsize(w);
    %calc PSNR Vs Lambda
    for n=1:length(Params.NN)
        curr=Params;                            
        curr.NN=curr.NN(n);                     values=struct2cell(curr);
        data_lambda(:,:,n,w)=StructPsnr(full_Data.ContextPsnr,fields,values,sigma_array,'Lambda');
        
    end
    %calc PSNR Vs NN
    for l=1:length(Params.Lambda)
        curr=Params;
        curr.Lambda=curr.Lambda(l);             values=struct2cell(curr);
        data_NN(:,:,l,w)=StructPsnr(full_Data.ContextPsnr,fields,values,sigma_array,'NN');
    end
    data_reference(w,1,:)=reshape(StructPsnr(full_Data.psnr,{'normalize','wsize'},values(4:5),sigma_array),1,1,[]);
end
data_lambda =reshape(data_lambda,l,length(sigma_array),[]); % dim1=lambda dim2=sigma dim3= NN X wsize
data_NN     =reshape(data_NN,n,length(sigma_array),[]);

ImproveLambda = shiftdim(data_lambda,2)-repelem(data_reference,n,l,1); %dim1=NN X wsize dim2=lambda %dim3=sigma
ImproveNN     = shiftdim(data_NN,2)-repelem(data_reference,l,n,1);

ImproveLambda(isnan(ImproveLambda))=0;
ImproveNN(isnan(ImproveNN))=0;

legend_lambda={};
for ind=1:(w*n)
    [nn,w1]=ind2sub([n,w],ind);
    legend_lambda={legend_lambda{:},strcat('NN:',num2str(Params.NN(nn)),'; W_1:',num2str(Arrays.wsize(w1)))};
end
legend_NN={};
for ind=1:(w*l)
    [lam,w1]=ind2sub([l,w],ind);
    legend_NN={legend_NN{:},strcat('\lambda:',num2str(Analysis.Arrays.Lambda(lam)),'; W_1:',num2str(Arrays.wsize(w1)))};
end
for f=ToPlot
    figure;
    subplot(2,1,1);
    plot(Analysis.Arrays.Lambda(1:l),ImproveLambda(:,:,f),'^-','MarkerSize',4); 
    title( strcat('\color{red} noise \sigma= ',num2str(sigma_array(f)),...
        '; CoOc Thr= ',num2str( Analysis.Arrays.CoOcThr(Thr_index)  ) )   )
    line(Analysis.Arrays.Lambda(1:l),0,'LineStyle','-.','Color','black')
    ylabel('Psnr improve Vs. \lambda');
    xlabel('\lambda')
    legend(legend_lambda,'Location','northeastoutside')
    subplot (2,1,2)
    plot(Analysis.Arrays.NN(1:n),ImproveNN(:,:,f),'*-','MarkerSize',4)
    line(Analysis.Arrays.NN(1:n),0,'LineStyle','-.','Color','black')

    ylabel('Psnr improve Vs. NN');
    xlabel('NN')
    legend(legend_NN,'Location','northeastoutside')
end


clearvars data_NN data_lambda curr