function [varargout]=PlotFromStruct(Method_1,Method_2,wsizes,varargin)
%% []=PlotData(psnr,RefPsnr,wsizes)
% takes struct psnr and plots it to graph

W=length(wsizes);
sigma=20:20:100;
if nargin>3
    handle=varargin{end};
    figure(handle)
    
else %create New plit
    handle=figure;
    for w=1:W
        ylabel(['window size=',num2str(wsizes(w))])
    end
xlabel('sigma');
subplot(W,1,1);legend('Method 1','Method 2')
ylabel('psnr')
end

for w=1:W

    [psnr_1]=struct2data (Method_1,'psnr.Lambda1.normalize0',wsizes(w));
    [psnr_2]=struct2data (Method_2,'psnr.Lambda1.normalize0',wsizes(w));
    
    subplot(W,1,w); hold on; grid on
    if isempty (psnr_2)
        plot(sigma,psnr_1,'--y*');
    elseif nargin>3
        plot(sigma,psnr_1,'--gd',sigma,psnr_2,'-.ms');
    else
        plot(sigma,psnr_1,'-b*',sigma,psnr_2,'-.r^');
    end
end
    
if nargout;varargout{1}=handle;end
end
%    Data_NN_5.psnr.Lambda1.normalize0.wsize7  
