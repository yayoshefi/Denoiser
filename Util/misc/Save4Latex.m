function []=Save4Latex(varargin)
% moves Images from the results directory 2 the latex figures

pnames = {'name'    'lambda'   'Thr'    'rule'   'CoOc',    'sigma'     'day'   'Method'    'loc'};
dflts =  {'',       0.1,          0,      4,      'MI',         50  ,    date()  'kmeans','C:\Users\yair-pc\Dropbox\Thesis code',};
[name,lambda,Thr,rule,CoOc,sigma,day,Method,location] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});

src=strcat(location,'\Results\',day,'\',name,'sigma',num2str(sigma),...
    '\rule #',num2str(rule),'_Context_comeans_',num2str(lambda,'%G'),...
    '_CoOc_',CoOc,'_Thr',num2str(Thr,'%G'),Method);
srcCoOc=strcat(location,'\Results\',day,'\',name,'sigma',num2str(sigma),'\CoOc.png');

trgt=strcat(location(1:end-5),'\Latex\figures\',day,'\',name,CoOc,'_sigma',num2str(sigma));

copyfile(src,trgt)
movefile(srcCoOc,trgt)


end