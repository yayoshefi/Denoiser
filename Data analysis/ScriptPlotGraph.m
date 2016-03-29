%% script for making Graphs

% change if necessary: struct name,strctfields,strctvalues
strctfields={'Lambda'   ,'NN'   ,'CoOcThr'  ,'normalize'    ,'wsize'};
strctvalues={ 1         ,7      ,1          ,0              ,9      }; %defualt stats

variable1='Lambda'; 
Var1_array=[1,2];
LogicVar1=ismember(strctfields,variable1);


data_context=[];    legendstring={};
for i=1:length(Var1_array)
strctvalues{LogicVar1}=Var1_array(i);

variable2='NN';         %can compare multiplie variables
Var2_array=[7,9];
LogicVar2=ismember(strctfields,variable2);
strctvalues{LogicVar2}=Var2_array;

[data_compare]=StructPsnr(full_Data.ContextPsnr,strctfields,strctvalues,sigma_array,variable2);
data_context=[data_context;data_compare];

%setting legend
    for m=Var2_array
        if variable1=='Lambda'
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

figure;plot (sigma_array,data_reference,'-.+k',sigma_array,data_context');
title('psnr rates for context clustering');legend('kmeans only',legendstring{:});
xlabel({cell2mat(strctfields);num2str(cell2mat(strctvalues))}); %FIX

gap=data_context-data_reference(ones(size(data_context,1),1),:);
figure;plot (sigma_array,gap,'-*');grid
title('diff psnr rates between context clustering and only visual clustering');legend(legendstring{:});

clearvars LogicVar1 LogicVar2 data_compare