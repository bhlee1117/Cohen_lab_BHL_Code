function sqrerror = fitgaussian2Dbgrnd_fixedwidth(parameters,width,data,X,Y)
%function sqrerror = fitgaussian2Dbgrnd_fixedwidth(parameters,width,data,X,Y)
%For use with fitminunc
%parameters(1) == mean(x)
%parameters(2) == mean(y)
%parameters(3) == amplitude of gaussian
%parameters(4) == background
%APF 2009-01-22

sqrerror = sum(sum((parameters(3)*(exp(-0.5*(X-parameters(1)).^2./(width^2)-0.5*(Y-parameters(2)).^2./(width^2))) + parameters(4) - data).^2));
