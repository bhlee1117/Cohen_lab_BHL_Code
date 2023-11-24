function [Y1, A, t, C, ci] = expplusconstfit(X, Y, A0, t0, C0, X1)
% 
% function [Y1, A, t, C, ci] = expplusconstfit(X, Y, A0, t0, C0, X1)
% fits an exponential plus a constant to the data Y(X).
% The exponential is given by A exp(-X/t)+ C
% Y1 is the fitted function, on the values specified by X1
% ci is the confidence interval in the fit
% 
% AEC 12/24/08; repaired AEC 2/28/2011

fun = inline('b(1)*exp(-X/b(2)) + b(3)', 'b', 'X');

[b R J] = nlinfit(X, Y, fun, [A0, t0, C0]);
ci = nlparci(b,R,J, 1e-12);
A = b(1);
t = b(2);
C = b(3);

Y1 = A * exp(-(X1/t))+ C;
