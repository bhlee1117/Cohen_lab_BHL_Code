% fitting to an exponential
function [Fcorrected,fit_trace]=PBcorr(Fmean)
nframes = length(Fmean);
n = (1:nframes)';

% the built-in library function for an exponential assumens the form 
% y = a*exp(b*x).  Values should be positive, so a > 0.  Photobleaching 
% will always be a decaying exponential, so b < 0.

% To get a good guess for the starting values, use a_start = max(y).
% To get a good gess at the exponential, make a linear approximation
% b ~ (y(end) - y(1)/max(a)).

Max = [Inf, 0];     % upper bound for [a,b]
Min = [0, -Inf];    % lower bound for [a,b]
Start = [max(Fmean), (Fmean(end) - Fmean(1))/max(Fmean)];   % starting guess for [a,b]

f = fit(n,Fmean,'exp1','StartPoint',Start,'Lower',Min,'Upper',Max);

fit_trace = f.a*exp(f.b*n);     % the fit function

% figure(2)
% plot(n,Fmean,'.',n,fit_trace,'r')
% legend('data','fit')

% correct for photobleaching
Fcorrected = Fmean./fit_trace*max(fit_trace);

% figure(3)
% subplot(2,1,1);plot(n,Fmean);legend('raw data');
% subplot(2,1,2);plot(n,Fcorrected);legend('corrected');


