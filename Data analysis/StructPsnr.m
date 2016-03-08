function [data]=StructPsnr(Strct,fields,values,sigma,varargin)
%% StructPsnr(Strct,fields,values,sigma)
% set fields as a cell array contaning the path is the structure
% values shoud be the parameter that are set in the path
% e.g. StructPsnr(Data_NN_5,{'Lambda','normalize','wsize'},{1,0,7},[10,20,30,40,60]);
%
% if we need to compare by one od the values, insert it as a cell array
% e.g data=StructPsnr(Data_NN_5,{'Lambda','normalize','wsize'},{1,0,[7,9]},[10,20,30,40,60],'wsize');

data=[];
comparefield=ismember(fields,varargin);
cmpvalues=values;

if sum(comparefield);arguments=length(values{comparefield});
else arguments=1;
end
for m=1:arguments
    if sum(comparefield);cmpvalues{comparefield}=values{comparefield}(m);
    end
    
%     pointer=Strct.psnr;
    pointer=Strct;
for idx=1:length(fields)
    substrct=strcat(fields{idx},num2str(cmpvalues{idx}));
    pointer=pointer.(substrct);
end
if isfield (pointer,strcat('sigma',num2str(sigma(1))))
    cell_data=struct2cell( pointer );
    if length (cell_data)>length(sigma);cell_data=cell_data(1:length(sigma));end
    data=[data;cell2mat(cell_data)'];
end

end
end