function [varargout]=CoOc_V1 (CoOc,visual,output,eps)
%% [varargout]=CoOc_V1 (CoOc,visual,output,eps)
% output = 'eps' \ 'entropy' \ 'both'

global Analysis
K=length(CoOc);

eps1=0.005;         eps2=0.02;          eps3=0.0001; 
if exist('eps','var');          eps1=eps;   end
epsNorm0=sum(CoOc(:)> 0 );      epsNorm1=sum(CoOc(:)>eps1); 
epsNorm2=sum(CoOc(:)>eps2);     epsNorm3=sum(CoOc(:)>eps3);

% Thresholding
% Thr=1.5*min(CoOc(CoOc~=0));
% CoOc(CoOc<Thr)=0;

% Entropy
LogCoOc=log2(CoOc);
LogCoOc(CoOc==0)=0;  % no nan resulting from inf*0;
H_row=-sum( CoOc.*LogCoOc,2 );
MatrixEntropy=mean(H_row);

if ~exist('output','var');output='none';end
switch lower(output)
    case 'eps'
        varargout={epsNorm0,epsNorm1,epsNorm2,epsNorm3};
    case 'entropy'
        varargout={MatrixEntropy};
    case 'both'
        varargout={MatrixEntropy,epsNorm0,epsNorm1,epsNorm2};
end

% Plot
if ~visual;    return
else
figure; colormap jet;
subplot(2,2,1); imagesc(log(CoOc+1));
title ('logaritmic Co-Occurence matrix');axis image
xlabel(['mean Entropy per row: ',num2str(MatrixEntropy)],'Color','red')

subplot(2,2,2);
sparsity=sum(CoOc>0,2);
plot(sparsity);title (strcat('using Thr: ','none'));       %num2str(Parameter.spatial.CoOcThr)));
grid on; ylabel('||CC_{k,l}||_{\epsilon}')
xlabel({'Labels K',...
    strcat('\color{blue} \epsilon=0; |Co-Oc|_{\0} : ',num2str(epsNorm0)),...
    strcat('\color{magenta} \epsilon=',num2str(eps1),'; |Co-Oc|_{\epsilon} : ',num2str(epsNorm1))});

subplot(2,2,[3,4]);
bar( 1:K ,CoOc(round(K/2),:) )
ylimit=min(1, 3*max(CoOc(:)) );
axis([1,K,0, ylimit])
title ( strcat('Co-Occurrence prob. for label ', num2str( round(K/2) )) );
if (isfield(Analysis,'K2')); line2=strcat( 'active context clusters ',num2str(Analysis.K2) );else line2='';end
xlabel ({strcat( 'amount of labels ',num2str(K) ),line2});axis ([1,K,0,1]); grid minor        
end

end