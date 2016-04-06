%% script for making Graphs to learn parameters

% change if necessary: struct name,strctfields,strctvalues
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

values=[data_reference',data_context'];
values=values-data_reference(ones(size(values,2),1),:)';

[data_reference]=StructPsnr(full_Data.psnr,{'normalize','wsize'},strctvalues(4:5),sigma_array);


figure;plot ([0,Var2_array],values,'-   *');
title('improvment Vs \lambda');%legend('kmeans only',legendstring{:});
xlabel(variable2); %FIX

% gap=data_context-data_reference(ones(size(data_context,1),1),:);
% figure;plot (sigma_array,gap);grid
% title('diff psnr rates between context clustering and only visual clustering');legend(legendstring{:});

clearvars LogicVar1 LogicVar2 data_compare