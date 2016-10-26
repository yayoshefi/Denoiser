function []=PrintDnoise (Output,result,varargin)

global Parameter Analysis
Method=Parameter.Method;
i=length(result);

%% De-noised images
h1=figure('name','Main Image after DeNoising');
subplot (1,2,1);
if verLessThan('matlab','8.4')
    % -- Code to run in MATLAB R2014a and earlier here --
    h1Num = num2str(h1);
else
    % -- Code to run in MATLAB R2014b and later here --
   h1Num = num2str(h1.Number);
end
if i==3
        ORACLE{1}=['  (',num2str(result{3}-result{1},'%1.3f'),')'];
        if result{2}-result{1}>0; sgn='+';else sgn='';end
        ORACLE{2}=['  (',sgn,num2str(result{2}-result{1},'%1.3f'),')'];
else    ORACLE{1}='';   ORACLE{2}='';
end
imshow(Output{1},[]); title(['psnr: ',num2str(result{1}),ORACLE{1}],'Color','r');
if Analysis.DebuggerMode && isfield(Analysis,'samp');DrawPixels (h1, Analysis.samp,true);end

xlabel({[Method,' main parameter ',num2str(Parameter.values.(Method))],...
    [' normalize parameter ',num2str(Parameter.normalize)]});

if i>1
ax_psnr2=subplot (1,2,2);
imshow(Output{2},[]); title(['psnr: ',num2str(result{2}),ORACLE{2}],'Color','b');
if Analysis.DebuggerMode && isfield(Analysis,'samp');DrawPixels (h1, Analysis.samp,true);end

xlabel({['Context: ',Parameter.Context,'. ',num2str( Analysis.K2),' clusters']})

end
dim = [0.01 0.65 0.3 0.3];
str=str2notation();
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
if Parameter.ORACLE
    annotation ('textbox',[0.4,0.85,0.1,0.1],'Color','b','FaceAlpha',0.2,...
        'String',strcat('ORACLE-',num2str(Analysis.ORACLE.level)));end

%% Labels images
figsize=[10 100 1100 600];
if nargin>2
    h2=figure('name','Lables','Position',figsize);
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
    subplot(1,2,1); 
    imshow(Lables1,[])
    xlabel(Method);
    if i>1
    subplot(1,2,2);
    imshow(Lables2,[])
    xlabel(['Context: ',Parameter.Context])
    xlabel(ax_psnr2,strcat('Context: ',Parameter.Context,num2str( length(unique(Lables2(:))) ),' clusters'))
    end
    colormap (Analysis.ColorMap);
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
if Parameter.ORACLE
    annotation ('textbox',[0.5,0.1,0.2,0.2],'Color','b','FaceAlpha',0.2...
        ,'String',strcat('ORACLE-',num2str(Analysis.ORACLE.level)));end
end
%% SAVING DATA
if Analysis.Save
    if isfield(Parameter,'ImageName');      prefix(Parameter.ImageName);
    else                                    prefix();
    end
    
    mkdir( strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix) );
    if ischar(result{1})
        saveas(h1,strcat(Parameter.location,'\Results\',date,'\',Parameter.ImageName,'_sigma',num2str(Parameter.sigma),...
            '\static\Inputs.png'));
    else
    saveas(h1,strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix,'\',h1Num,Analysis.ShortPrefix,'_D-noisedImages.png'))
    end
    if exist('h2','var')         % print Lables if there is
        saveas(h2,strcat(Parameter.location,'\Results\',date,'\',Analysis.DirSuffix,'\',h2Num,Analysis.ShortPrefix,'_Lables.png'))
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

function str=str2notation()
global Parameter
switch Parameter.CoOc.Type
    case 'MI';      Type='Mutual info';
    case 'PMI';     Type='PointWise Mutual info';
    case 'CC';      Type='Conditional';
    case 'JP';      Type='Joint Prob.';
end
str={['\bf Context: ',Parameter.Context],...
    ['\rm \lambda: ',num2str(Parameter.spatial.lambda,'%G'),' NN:',num2str(Parameter.spatial.NN)]...
    ,[Type,' Shrinkage  %=',num2str(Parameter.CoOc.ShrinkPer)],' Thr@',num2str(Parameter.CoOc.Thr,'%2.1G'),...
    ['W_1: ',num2str(Parameter.wsize2^0.5),'   \sigma: ',num2str(Parameter.sigma)],...
    ['Update rule: ',num2str(Parameter.spatial.UpdateRule),' with ',Parameter.CoOc.AssginType,' assignment'],...
    'reserved..'};

end

function []= prefix(name)
global Parameter Analysis
% in PrintDnoise there is a subfolder for diffrent noise
if ~exist('name','var');name='';end
Analysis.DirSuffix=strcat(name,'_sigma',num2str(Parameter.sigma),'\',Parameter.description,...
    '_Context_',Parameter.Context,...
    '_',num2str(Parameter.spatial.lambda,'%G'),...
    '_CoOc_',Parameter.CoOc.Type,'_Thr',num2str(Parameter.CoOc.Thr,'%G'),...
    Parameter.Method);
    %,'_norm',num2str(Parameter.normalize),'_metric_',Parameter.metric);

Analysis.ShortPrefix=strcat('Rule-',num2str(Parameter.spatial.UpdateRule),Parameter.CoOc.AssginType,...
    '_W1-',num2str(sqrt(Parameter.wsize2)),'_Context-',Parameter.Context,...
    '-',num2str(Parameter.spatial.lambda,'%G'),'_NN-',num2str(Parameter.spatial.NN),';');
end


