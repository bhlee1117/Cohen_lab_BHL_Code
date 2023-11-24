function [y_fit t_consts coeffY]  = expfitDM_2(x,y,xs,t_guess)
% finds a least-squares fit of n exponentials (plus a const)
% xs is x points at whic to sample for output (allowing extrapolation etc)
% x, y, xs should be column vectors.  t_guess should be a row vector.
% 4-19-2012 Now making it fit multiple curves (all to the same
% exponentials)

xs = xs-x(1);
x  = x-x(1); % Doing this to avoid the silly size of exponentials if there's a time offset

errfun = @(param) mean(mean((y - exp(- x*[0, param])*fit_curves(y, exp(-x*[0, param]))).^2  , 1),2);
t_consts = 1./fminsearch(errfun,1./t_guess);

coeffY = fit_curves(y,exp(-x*[0, 1./t_consts]));
y_fit =   exp(-xs*[0, 1./t_consts])*coeffY;

end