function [data]=struct2data (strct,field,windowsize)
%% data is a row vector taken from: name.field.wsize'windowsize'...

if isstruct (strct)
    cell_data=struct2cell( strct.(field).(['wsize',num2str(windowsize)]) );
    data=cell2mat(cell_data)';
else
    data=[];
% cell_centers=struct2cell( name.centers.(['wsize',num2str(wsizes(w))]) );
% centersCount=cell2mat(cell_centers)';
end