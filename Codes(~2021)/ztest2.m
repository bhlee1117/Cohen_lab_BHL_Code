function zval = ztest2(x,y)
% x - sample 1
% y - sample 2
% mux - mean of x to be used in difference. If you are not testing a
% specific difference, use mux = 0
% muy - mean of y. If you are not testing a
% specific difference, use muy = 0
% varx - variance of x
% vary - variance of y
    Nx = length(x);
    Ny = length(y);
    p=(sum(x)+sum(y))/(Nx+Ny);
  zval = (mean(x)-mean(y))/sqrt(p*(1-p)*(1/Nx+1/Ny));
  end