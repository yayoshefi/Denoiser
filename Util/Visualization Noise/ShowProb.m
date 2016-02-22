function []=ShowProb (affinity,samp)
%% []=ShowProb (affinity,samp)
% the function plots the PDF distribution of a sample of the data based on
% the given full affinity.
% affinity is [n By K]-n data points over K options
%
% case samp is a scalar, the the plot for samp random data points
% case samp is a vector, the plot is for the given data points in samp


[n,K]=size(affinity);
if isscalar (samp) ; samp=randperm(n,samp);end

data=affinity(samp,:);
figure;
bar(data); colormap jet
title (strcat ('prob. for pixels: ',num2str(samp)));
c=colorbar;c.Label.String='Labels';
xlabel ('pixel');ylabel('Probabilty')

% 
% colors={'b','r','k','g','m','y'};
% subplot(1,2,2);b=bar(1:K,data','histc');title (strcat ('prob. for pixels: ',num2str(samp)));
% xlabel ('Label');ylabel('Probabilty')
% for i=1:length(samp)
%     b(i).FaceColor=colors{i}; b(i).EdgeColor=colors{i};
% end
end