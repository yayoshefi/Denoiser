function I=tf_idf (H)

[n,K]=size (H);

tf=H;
idf=sum(logical(H),2);
idf=log(1+idf);

I=tf.*idf(:,ones(1,K));

end

%I_pixel=sum(I,2); 

%how can in meassure the information in H(i) regarding all the data and not
%just one of the labels?