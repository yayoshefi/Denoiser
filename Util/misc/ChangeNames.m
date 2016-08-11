%% This Script changes the name in the path 

PthDir='C:\Users\Yair\Dropbox\Thesis\Latex\figures\11-Aug-2016\texture_CC_sigma50';
NewDir=PthDir;

ind=strfind(NewDir,'-');
for i=length(ind):-1:1
    NewDir(ind(i))='';
end
ind=strfind(NewDir,'_');
for i=length(ind):-1:1
    NewDir(ind(i))='';
end
mkdir(NewDir);

fils=dir(PthDir);
for f=1:length(fils)
    if ( strfind(fils(f).name,'Params') )
        delete ([PthDir,'\',fils(f).name])
    end
end
%{
for f=1:length(fils)
    for r=1:7
        if ( strfind(fils(f).name,['Rule-',num2str(r)]) )
            if (strfind(fils(f).name,'Lables'))
                NewName=strcat('Rule',num2str(r),'Labels');
                movefile([PthDir,'\',fils(f).name],[NewDir,'\',NewName,'.png'])
            elseif (strfind(fils(f).name,'noised'))
                NewName=strcat('Rule',num2str(r),'denoised');
                movefile([PthDir,'\',fils(f).name],[NewDir,'\',NewName,'.png'])
            end
            
        end
    end
    
end
%}
% movefile([PthDir,'\CoOc.png'],[NewDir,'\CoOc.png']);
movefile(PthDir,NewDir);

