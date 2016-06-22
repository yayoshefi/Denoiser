function []=ShowProb (affinity,samp)
%% []=ShowProb (affinity,samp)
% the function plots the PDF distribution of a sample of the data based on
% the given full affinity.
% affinity is [n By K]-n data points over K options
%
% to compare more than one affinity, use the 3rd dim
%
% case samp is a scalar, the the plot for samp random data points
% case samp is a vector, the plot is for the given data points in samp

[n,K,P]=size(affinity);
if isscalar (samp) ; samp=randperm(n,samp);end
clrs={'red','green','blue','cyan','magenta'};
str=[];
data=affinity(samp,:,:);
figure;
for p=1:P
    subplot(P,1,p);
    bar(data(:,:,p)); colormap jet;grid on;
    c=colorbar;c.Label.String='Labels';
    str=strcat(str,strcat(' \color{',clrs{mod(p-1,5)+1},'}',num2str(samp(p)) ));
%     if p==1;title (strcat ('prob. for pixels: ',num2str(samp)));        end    
    if p==P;xlabel ('pixel');ylabel('Probabilty');                      end
end
subplot(P,1,1);
title (strcat ('prob. for pixels: ',str ));
% 
% colors={'b','r','k','g','m','y'};
% subplot(1,2,2);b=bar(1:K,data','histc');title (strcat ('prob. for pixels: ',num2str(samp)));
% xlabel ('Label');ylabel('Probabilty')
% for i=1:length(samp)
%     b(i).FaceColor=colors{i}; b(i).EdgeColor=colors{i};
% end
end