function sqrerror = fitgaussian2Dbgrnd(parameters,data,X,Y)
%function sqrerror = fitgaussian2Dbgrnd(parameters,data,X,Y)
%For use with fminunc
%parameters(1) == mean(x)
%parameters(2) == mean(y)
%parameters(3) == width (x and y)
%parameters(4) == amplitude of gaussian
%parameters(5) == background
%APF 2009-01-22

sqrerror = sum(sum((parameters(4)*(exp(-0.5*(X-parameters(1)).^2./(parameters(3)^2)-0.5*(Y-parameters(2)).^2./(parameters(3)^2))) + parameters(5) - data).^2));
