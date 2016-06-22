function []=DrawPixels (h, pxl,space)
%% DrawPixels (h, pxl)
% Input:
%   h-      handle to image
%   pxl-    index pixels in a cropped image (after imfilter)
%   space-  boolean field. true-visual space; false- label space
global Parameter
clrs={'r','g','b','c','m'};

M=(sqrt(Parameter.wsize2)-1)/2;
[i,j]=ind2sub([Parameter.row-M,Parameter.col-M],pxl);
if space; i=i+M;j=j+M;  end        % change to large image coordinates

r=3;                %radius of the circle
theta = 0:pi/50:2*pi;
x = r * ( ones(size(pxl))' * cos(theta) ) + j'                   * ones (size(theta));
y = r * ( ones(size(pxl))' * sin(theta) ) + (Parameter.row - i)' * ones (size(theta));

figure(h); hold on
for p=1:length(pxl)
    plot(x(p,:)', y(p,:)',clrs{p},'LineWidth',0.9);
    plot(j(p), Parameter.row-i(p),strcat('+',clrs{p}),'MarkerSize',2,'LineStyle','none');
end

end