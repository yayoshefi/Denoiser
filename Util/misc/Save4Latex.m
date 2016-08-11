function []=Save4Latex(varargin)

pnames = {'name'    'lambda'   'Thr'    'rule'   'CoOc',    'sigma'     'day'};
dflts =  {'',       0.1,          0,      4,      'MI',         50  ,    date()};
[name,lambda,Thr,rule,CoOc,sigma,day] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});

src=strcat('C:\Users\Yair\Dropbox\Thesis code\Results\',day,'\',name,...
    'sigma',num2str(sigma),'\rule #',num2str(rule),'_Context_comeans_',num2str(lambda,'%G'),...
    '_CoOc_',CoOc,'_Thr',num2str(Thr,'%G'),'kmeans');
srcCoOc=strcat('C:\Users\Yair\Dropbox\Thesis code\Results\',day,'\',name,...
    'sigma',num2str(sigma),'\CoOc.png');

trgt=strcat('C:\Users\Yair\Dropbox\Thesis\Latex\figures\',day,'\',name,CoOc,'_sigma',num2str(sigma));

copyfile(src,trgt)
movefile(srcCoOc,trgt)


end