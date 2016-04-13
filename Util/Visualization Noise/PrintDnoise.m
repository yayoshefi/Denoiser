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
% ylabel (['# ',num2str(Analysis.K),' clusters'])
xlabel({[Method,' main parameter ',num2str(Parameter.values.(Method))],...
    [' normalize parameter ',num2str(Parameter.normalize)]});

if i>1
ax_psnr2=subplot (1,i,i);
imshow(Output{2},[]); title(['psnr: ',num2str(result{2})],'Color','b');
% ylabel (['Context, using  ',num2str(Parameter.Spectral.clustrsNUM),' clusters'])
xlabel({['Context: ',Parameter.Context,num2str(Parameter.Spectral.clustrsNUM),' clusters']})
%     ['Method used: ',Parameter.Context,' Context lambda: ',num2str(Parameter.Spatil.lambda,'%G')]...
%     ,['windows size is: ',num2str(Parameter.wsize2^0.5),' noise std is: ',num2str(Parameter.sigma) ]})%,...
%     ['lambda: ',num2str(lambda)]});
end
dim = [0.01 0.65 0.3 0.3];
str={['\bf Context: ',Parameter.Context],...
    ['\rm \lambda: ',num2str(Parameter.Spatil.lambda,'%G'),' NN:',num2str(Parameter.Spatil.NN)]...
    ,['|CC|_{\epsilon}   \epsilon=',num2str(Parameter.Spatil.CoOcThr,'%G')],...
    ['W_1: ',num2str(Parameter.wsize2^0.5),'   \sigma: ',num2str(Parameter.sigma)],...
    'reserved..'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');


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
    subplot(1,i,1); 
    imshow(Lables1,[])
    xlabel(Method);
    if i>1
    subplot(1,i,2);
    imshow(Lables2,[])
    xlabel(['Context: ',Parameter.Context])
    xlabel(ax_psnr2,strcat('Context: ',Parameter.Context,num2str( length(unique(Lables2(:))) ),' clusters'))
    end
    colormap (Analysis.ColorMap);
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
end

if Analysis.Save
    prefix();
    mkdir( strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix) );
    
    saveas(h1,strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix,'\',h1Num,Analysis.ShortPrefix,'_D-noisedImages.jpg'))
    if exist('h2','var')         % print Lables if there is
        saveas(h2,strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix,'\',h2Num,Analysis.ShortPrefix,'_Lables.jpg'))
    end

    % write parameters to text file
    fileID = fopen(strcat(Parameter.location,'\Results\',date(),'\',Analysis.DirSuffix,'\','Params_',Analysis.ShortPrefix,'.txt'),'w');
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

function []= prefix()
global Parameter Analysis
% in PrintDnoise there is a subfolder for diffrent noise
Analysis.DirSuffix=strcat('sigma',num2str(Parameter.sigma),'\',Parameter.description,'Window',num2str(sqrt(Parameter.wsize2)),...
    '_Context_',Parameter.Context,...
    '_',num2str(Parameter.Spatil.lambda,'%G'),'NN_',num2str(Parameter.Spatil.NN),...
    'CoOcThr',num2str(Parameter.Spatil.CoOcThr,'%G'),...
    Parameter.Method,...
    '_norm',num2str(Parameter.normalize),'_metric_',Parameter.metric);

Analysis.ShortPrefix=strcat('W1-',num2str(sqrt(Parameter.wsize2)),'_Context-',Parameter.Context,...
    '-',num2str(Parameter.Spatil.lambda,'%G'),'_NN-',num2str(Parameter.Spatil.NN),';');
end


