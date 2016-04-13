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
Params.CoOcThr=Params.CoOcThr(1);   Params.normalize=Params.normalize(1);


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
    plot(Analysis.Arrays.Lambda(1:l),ImproveLambda(:,:,f),'^-','MarkerSize',4); title(strcat('\color{red} sigma=',num2str(sigma_array(f))))
    line(Analysis.Arrays.Lambda(1:l),0,'LineStyle','-.','Color','black')
    ylabel('Psnr improve Vs. \lambda');
    xlabel('\lambda')
    legend(legend_lambda)
    subplot (2,1,2)
    plot(Analysis.Arrays.NN(1:n),ImproveNN(:,:,f),'*-','MarkerSize',4)
    line(Analysis.Arrays.NN(1:n),0,'LineStyle','-.','Color','black')

    ylabel('Psnr improve Vs. NN');
    xlabel('NN')
    legend(legend_NN)
end


%{
strctfields={'Lambda'   ,'NN'   ,'CoOcThr'  ,'normalize'    ,'wsize'};
strctvalues={ 1         ,3      ,2          ,0              ,5      }; %defualt stats

variable1='NN'; 
Var1_array=[3];
LogicVar1=ismember(strctfields,variable1);


data_context=[];    legendstring={};
for i=1:length(Var1_array)
strctvalues{LogicVar1}=Var1_array(i);

variable2='Lambda';         %can compare multiplie variables
Var2_array=[1,2,3];
LogicVar2=ismember(strctfields,variable2);
strctvalues{LogicVar2}=Var2_array;

[data_compare]=StructPsnr(full_Data.ContextPsnr,strctfields,strctvalues,sigma_array,variable2);
data_context=[data_context;data_compare];

%setting legend
    for m=Var2_array
        if strcmp(variable1,'Lambda')
            legendstring={legendstring{:},strcat(variable1,'='...
                ,num2str(lambda_array(Var1_array(i))),';',variable2,'=',num2str(m))};
        elseif variable2=='Lambda'
            legendstring={legendstring{:},strcat(variable1,'='...
                ,num2str(Var1_array(i)),';',variable2,'=',num2str(lambda_array(m)))};
        else
            legendstring={legendstring{:},strcat(variable1,'='...
                ,num2str(Var1_array(i)),';',variable2,'=',num2str(m))};
        end

    end
end
[data_reference]=StructPsnr(full_Data.psnr,{'normalize','wsize'},strctvalues(4:5),sigma_array);

values=[data_reference',data_context'];
values=values-data_reference(ones(size(values,2),1),:)';


figure;plot ([0,Var2_array],values,'-   *');
title('improvment Vs \lambda');%legend('kmeans only',legendstring{:});
xlabel(variable2); %FIX

% gap=data_context-data_reference(ones(size(data_context,1),1),:);
% figure;plot (sigma_array,gap);grid
% title('diff psnr rates between context clustering and only visual clustering');legend(legendstring{:});
%}
clearvars data_NN data_lambda curr