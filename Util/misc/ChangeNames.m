%% This Script changes the name in the path 

PthDir='C:\Users\Yair\Dropbox\Thesis\Latex\figures\24-Aug-2016\chateau_CC_sigma50';
NewDir=PthDir;

NewDir=strrep(NewDir,'_','');
NewDir=strrep(NewDir,'-','');
mkdir(NewDir);

% folders=dir(PthDir);
% if folders(end).isdir;    cnt_folders=length(folders);  %how to find if the is one image or folders?
% else                    cnt_folders=1;              end
% 
% for fold=3:cnt_folders
%     if folders(fold).isdir; files=dir([PthDir,'\',folders(fold).name]);
%     else files=dir(PthDir);                                             end

files=dir(PthDir);    

for f=1:length(files)
    for r=1:7
        if ( strfind(files(f).name,['Rule-',num2str(r)]) )
            if (strfind(files(f).name,'Lables'))
                NewName=strcat('Rule',num2str(r),'Labels');
                movefile([PthDir,'\',files(f).name],[NewDir,'\',NewName,'.png'])
            elseif (strfind(files(f).name,'noised'))
                NewName=strcat('Rule',num2str(r),'denoised');
                movefile([PthDir,'\',files(f).name],[NewDir,'\',NewName,'.png'])
            end
            
        end
    end
    
end
% end
%movefile(PthDir,NewDir);

static={'CoOc','ORACLE','ORACLEdenosing','ORACLECoOc','DM3Ddenosing','Inputs'};
for s=static
    movefile([PthDir,'\',s{:},'.png'],[NewDir,'\',s{:},'.png']);
end


