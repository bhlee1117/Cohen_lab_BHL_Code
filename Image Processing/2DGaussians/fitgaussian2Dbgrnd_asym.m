function sqrerror = fitgaussian2Dbgrnd_asym(parameters,data,X,Y)
%function sqrerror = fitgaussian2Dbgrnd_asym(parameters,data,X,Y)
%For use with fminunc
%parameters(1) == mean(x)
%parameters(2) == mean(y)
%parameters(3) == width(x)
%parameters(4) == width(y)
%parameters(5) == amplitude of gaussian
%parameters(6) == background
%APF 2009-11-17

sqrerror = sum(sum((parameters(5)*(exp(-0.5*(X-parameters(1)).^2./(parameters(3)^2)-0.5*(Y-parameters(2)).^2./(parameters(4)^2))) + parameters(6) - data).^2));
