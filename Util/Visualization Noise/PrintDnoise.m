function []=PrintDnoise (Output,result,varargin)

global Parameter Analysis
Method=Parameter.Method;
i=length(result);


h1=figure('name','Main Image after DeNoising');
subplot (1,i,1);
if verLessThan('matlab','8.4')
    % -- Code to run in MATLAB R2014a and earlier here --
    h1Num = num2str(h1);
else
    % -- Code to run in MATLAB R2014b and later here --
   h1Num = num2str(h1.Number);
end

imshow(Output{1},[]); title(['psnr: ',num2str(result{1})],'Color','r');
ylabel (['# ',num2str(Analysis.K),' clusters'])
xlabel({['Method used: ',Method,' main parameter ',num2str(Parameter.values.(Method)),...
    ' normalize parameter ',num2str(Parameter.normalize)],...
    ['windows size is: ',num2str(Parameter.wsize2^0.5),' noise std is: ',...
    num2str(Parameter.sigma)]});

if i>1
subplot (1,i,i);
imshow(Output{2},[]); title(['psnr: ',num2str(result{2})],'Color','b');
ylabel (['Context, using  ',num2str(Parameter.Spectral.clustrsNUM),' clusters'])
xlabel({['Method used: ',Parameter.Context,' Context lambda: ',num2str(Parameter.Spatil.lambda,'%G')]...
    ,['windows size is: ',num2str(Parameter.wsize2^0.5),' noise std is: ',num2str(Parameter.sigma) ]})%,...
%     ['lambda: ',num2str(lambda)]});
end



if nargin>2
    h2=figure('name','Lables');
    if verLessThan('matlab','8.4')
        % -- Code to run in MATLAB R2014a and earlier here --
        h2Num = num2str(h2);
    else
        % -- Code to run in MATLAB R2014b and later here --
       h2Num = num2str(h2.Number);
    end    
    wsize=Parameter.wsize2^0.5;
    [r1,c1]=size(Output{1});
    [r2,c2]=size(Output{2});    
    
    Lables1=col2im(varargin{1},[wsize,wsize],[r1,c1],'sliding');
    Lables2=col2im(varargin{2},[wsize,wsize],[r2,c2],'sliding');
    subplot(1,i,1); colormap jet
    imagesc(Lables1)
    if i>1
    subplot(1,i,2);
    imagesc(Lables2)
    end

end

if Analysis.Save
    str=strcat(Parameter.description,Parameter.Method,'_metric_',...
        Parameter.metric,'_Context_',Parameter.Context,...
    '_',num2str(Parameter.Spatil.lambda,'%G'),'_normalize',num2str(Parameter.normalize),...
    '_sigma',num2str(Parameter.sigma));
    
    mkdir( strcat(Parameter.location,'\Results\',date,'\',str) );
    
    saveas(h1,strcat('Results\',date,'\',str,'\',h1Num,'_D-noisedImages.jpg'))
    if exist('h2','var')         % print Lables if there is
        saveas(h2,strcat('Results\',date,'\',str,'\',h2Num,'_Lables.jpg'))
    end

    % write parameters to text file
    fileID = fopen(strcat('Results\',date(),'\',str,'\','Parameter.txt'),'w');
    fprintf(fileID,ParameterText(Parameter));
    fclose(fileID);
end
end

function txt=ParameterText(strct)

% txt=strcat('Method used :',Parameter.Method);
txt=strcat(inputname(1),':  ',date(),'\r\n');
values=struct2cell(strct);
names = fieldnames(strct);
for i=1:length(values)
    if isstruct(values{i})
        txt=strcat( txt,'\r\n',names{i},ParameterText( strct.(names{i}) ),'\r\n' );
    else
        if ~strcmp(names{i},'B')
            txt=strcat(txt,num2str(i),') ',num2str(names{i}),':  ',num2str(values{i}),'\r\n');end
%             num2str(values{i}),'\r\n');
        
    end
end
end


