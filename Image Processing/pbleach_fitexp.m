function [fitted_curves, fitted_params, intensity_corrected] = ...
    pbleach_fitexp(intensity, Nskip, Nstop)
% fitted_curves = pbleach_fitexp(intensity)
% Fits an exponential decay to data in each column. Returns the fit curve.
% fittedcurves = pbleach_fitexp(intensity, Nskip) Skips Nskip initial
% samples from consideration to fit the curve.
% fitted_curves = pbleach_fitexp(intensity, Nskip, Nstop) Skips Nskip
% initial samples, and ignores samples after Nstop from consideration to
% fit the curve. 
% [fitted_curves, fitted_params] = pbleach_fitexp(...) Also returns the
% fitted parameters in the format [offs; amps; taus] where each column of
% intensity has been fitted with an exponential of the form
%     amp * exp(-t/tau) + off
% Each column of fitted_params has the parameters fitted for each column of
% intensity. 
% [fitted_curves, fitted_params, intensity_corrected] = pbleach_fitexp(...)
% Also returns intensity corrected traces in each column, after dividing
% input data by the fitted curves, and restoring the mean of each column.  
%
%     Example use: 
% 
%     rng(2019)
%     data = 50*(exp(-(1:500)'./(100*rand(1,3))).*randn(1,3) + rand(1,3)) + randn(500,3); 
%     fitdata = pbleach_fitexp(data);
%     figure
%     plot(data)
%     hold on
%     plot(fitdata)
%
% 
% 2019 Vicente Parot
% Cohen Lab - Harvard University
% 
    [Nsamples, Ntraces, ~] = size(intensity);
    if ~exist('Nstop','var'), Nstop = Nsamples; end
    if ~exist('Nskip','var'), Nskip = 0; end
    t = (0:Nsamples-1)';
    q = intensity(1+Nskip:Nstop,:);
    offs = mean(q(ceil(.75*end):end,:));
    amps = mean(q(1:ceil(.25*end),:))-offs;
    taus = Nsamples.*ones(1,Ntraces)/5;
    for it = 1:Ntraces
        dintens = double(intensity(Nskip+1:Nstop,it));
        expf = @(t,v) v(1,:) + (v(2,:)).*(exp(-t./v(3,:)));
        objf = @(v) sum(sum((dintens - expf(t(Nskip+1:Nstop),v)).^2));
        initialparams = [offs; amps; taus];
        fitted_params(:,it) = fminsearch(objf,initialparams(:,it));
        fitted_curves(:,it) = expf(t,fitted_params(:,it));
    end
    if nargout > 1
        intensity_corrected = intensity./fitted_curves.*mean(intensity);
    end
