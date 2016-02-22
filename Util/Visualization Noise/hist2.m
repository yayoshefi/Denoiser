function [h,xtics,ytics]=hist2(x,y,xbins,ybins,xlhb,ylhb)
%function [h,xtics,ytics]=hist2(x,y,xbins,ybins,xlhb,ylhb)
% 
%   generates 2D histogram of x an y
%  xbins,ybins are the number of bins in the histogram
%  xlhb=[lb,hb] and ylhb=[lb,hb] (optional) specify the lb and hb of each variable 

if (length(x)~=length(y))
    error('length(x) should be equal to length(y)');
end;

bins=xbins;
lb=min([x(:) ; y(:)] ); hb=max([x(:) ; y(:)]);
if exist('xlhb'),
    lb=xlhb(1); hb=xlhb(2);
end;
t1=[lb:(hb-lb)/bins:hb]';     
factor1=bins/(hb-lb);
shift1=-factor1*lb+1;

% bins=ybins;
% lb=min(y(:)); hb=max(y(:));
% if exist('ylhb'),
%     lb=ylhb(1); hb=ylhb(2);
% end;
% t2=[lb:(hb-lb)/bins:hb]';     
% factor2=bins/(hb-lb);
% shift2=-factor2*lb+1;

t2=t1;
factor2=factor1;
shift2=shift1;


h=zeros(xbins+1,ybins+1);
x=x(:);
y=y(:);

for i=1:size(x,1),
            s1=round(x(i)*factor1+shift1);
            s2=round(y(i)*factor2+shift2);
            
            h(s1,s2)=h(s1,s2)+1;
            
     end;
%  end;
 
 xtics=t1;
 ytics=t2;
 
 