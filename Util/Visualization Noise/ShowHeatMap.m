function [HeatMap] = ShowHeatMap (Original, Denoised)
%% [HeatMap] = ShowHeatMap (Original, Denoised)
% the functio shows the diffrence between the Original image and a denoised
% version. one can learn the strengths of the denoised algorithm from such
% compare.

HeatMap = abs( im2double(Original) - im2double(Denoised) );

if isvector(HeatMap)
    L=sqrt(length(HeatMap));
    HeatMap=reshape(HeatMap,L,L);
end

figure;
imshow(HeatMap,[]); colormap jet;