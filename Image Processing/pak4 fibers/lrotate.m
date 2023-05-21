function[xf, yf]=lrotate(x,y, a, x0, y0)
% x and Y are the input vectors whose line graph to be rotated.
% a: angle of rotation
%x0, y0 the x and y coordinate of the point about which the line is to be rotated
% returns: xf and yf which forms the rotated graph
%
n=length(x);
xf=NaN(size(x));
yf=NaN(size(y));
for i=1:n
 x1=x(i);
 y1=y(i);
 r=sqrt((y1-y0)^2+(x1-x0)^2);
 if((x1>=x0)&&(y1>=y0))
     b=atan2d((y1-y0),(x1-x0));
     x2=x0+r*cosd(a+b);
     y2=y0+r*sind(a+b);
 elseif((x1<x0)&&(y1>=y0))
     b=atan2d((y1-y0),(x0-x1));
     x2=x0-r*cosd(b-a);
     y2=y0+r*sind(b-a);
 elseif ((x1<x0)&&(y1<y0))
     b=atan2d((y0-y1),(x0-x1));
     x2=x0-r*cosd(b+a);
     y2=y0-r*sind(b+a);
 elseif ((x1>=x0)&&(y1<y0))
     b=atan2d((y0-y1),(x1-x0));
     x2=x0+r*cosd(b-a);
     y2=y0-r*sind(b-a);
 end
 xf(i)=x2;
 yf(i)=y2;
end
end
