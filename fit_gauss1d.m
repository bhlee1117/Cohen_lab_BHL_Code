function [pfit yhat R2]=fit_gauss1d(x,y)
% p = [A, mu, sigma, y0]

gauss = @(p,x) p(1)*exp(-((x-p(2)).^2)./(2*p(3)^2)) + p(4);
% p = [A, mu, sigma, y0]

y0_init = median(y);              % robust background
A_init  = max(y) - y0_init;
[~,imax] = max(y);
mu_init = x(imax);
sigma_init = (max(x)-min(x))/6;   % fallback guess

p0 = [A_init, mu_init, sigma_init, y0_init];

% --- bounds: sigma>0
lb = [-Inf, min(x), 1e-9, -Inf];   % allow negative amplitude/background if needed; adjust if you want positivity
ub = [ Inf, max(x), Inf,  Inf];

% --- perform fit
opts = optimoptions('lsqcurvefit','Display','off','MaxFunctionEvaluations',1e4);
[pfit,resnorm,~,exitflag,output,lambda,jac] = lsqcurvefit(gauss,p0,x,y,lb,ub,opts);

A = pfit(1); mu = pfit(2); sigma = pfit(3); y0 = pfit(4);
FWHM = 2*sqrt(2*log(2))*sigma;

% --- diagnostics & plot
yhat = gauss(pfit,x);
resid = y - yhat;
SSR = sum(resid.^2);
SST = sum((y - mean(y)).^2);
R2 = 1 - SSR/SST;
