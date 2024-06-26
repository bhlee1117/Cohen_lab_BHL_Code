function [z] = fitgaussian1D(p,v,x);

%cx = p(1);
%wx = p(2);
%amp = p(3);

zx = p(3)*exp(-0.5*(x-p(1)).^2./(p(2)^2)) - v;

z = sum(zx.^2);
