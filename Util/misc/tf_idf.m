function [I, PixelTotalI]=tf_idf (H)

[n,K]=size (H);

tf=H;
idf=sum(logical(H),1);
idf=log(1+n./idf);

I=tf.*idf(ones(n,1),:);

PixelTotalI=sum(I,2);
end

%I_pixel=sum(I,2); 

%how can in meassure the information in H(i) regarding all the data and not
%just one of the labels?