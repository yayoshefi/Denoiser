function []=DrawPixels (h, pxl,space)
%% DrawPixels (h, pxl)
% Input:
%   h-      handle to image
%   pxl-    index pixels in a cropped image (after imf2col) 
%           OR      [x1 ,y1 ;
%                    x2 ,y2  ] values
%   each pixel to print is at the 1st dim
%   in the Original size Image
%   space-  boolean field. true-visual space; false- label space
global Parameter
clrs={'r','g','b','c','m'};
M=(sqrt(Parameter.wsize2)-1)/2;

if size(pxl,1)>1 % pixel given in [X,Y] pairs from the Orignal Image, NEEDE at least 2 sets
    i=pxl(:,1)-M;   j=pxl(:,2)-M;
else
    pxl = pxl(:); % make sure you have a column vector
    [i,j]=ind2sub([Parameter.row-2*M,Parameter.col-2*M],pxl);
end

if space; i=i+M;j=j+M;  end        % change to large image coordinates

r=3;                %radius of the circle
theta = 0:pi/50:2*pi;
x = r * ( ones(size(pxl,1),1) * cos(theta) ) + ( j) * ones (size(theta));
y = r * ( ones(size(pxl,1),1) * sin(theta) ) + ( i) * ones (size(theta));

figure(h); hold on
for p=1:size(pxl,1)
    plot(x(p,:)', y(p,:)',clrs{p},'LineWidth',0.9);
    plot(j(p), i(p),strcat('+',clrs{p}),'MarkerSize',2,'LineStyle','none');
end

end