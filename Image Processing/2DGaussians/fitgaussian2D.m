function [z] = fitgaussian2D(p,m,X,Y);

%cx = p(1);
%cy = p(2);
%wx = p(3);
%wy = p(4);
%amp = p(5);

ztmp = p(5)*(exp(-0.5*(X-p(1)).^2./(p(3)^2)-0.5*(Y-p(2)).^2./(p(4)^2))) - m;

z = sum(sum(ztmp.^2));
