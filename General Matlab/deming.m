function [b0, b1] = deming(x, y, delta)
% function [b0, b1] = deming(x, y, delta)
% Calculate the Deming regression of y and x, with delta the ratio of the expected squared errors (y versus x)
% Minimizes the total error in the regression y = b0 + b1*x

s = cov(x,y);
b1 = (s(2,2) - delta*s(1,1) + sqrt((s(2,2) - delta*s(1,1))^2 + 4*delta*(s(1,2)^2)))/(2*s(1,2));
b0 = mean(y) - b1*mean(x);
