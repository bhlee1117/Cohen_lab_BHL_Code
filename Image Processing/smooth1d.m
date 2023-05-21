function out = smooth1d(in,sigma,boundary)
% gaussian lowpass smoothing filter with default sigma 15
%
% 2016 Vicente Parot
% Cohen Lab - Harvard University
%
if ~exist('sigma','var')
    sigma = 15;
end
if ~exist('boundary','var')
    boundary = 'none';
end
switch boundary
    case 'none'
        ns = size(in,1);
        gw = gausswin(floor(ns/2)*2+1,sigma);
        gw = gw(1:ns);
        out = ifft(bsxfun(@times,fft(in),ifftshift(gw)));
    case 'replicate'
        pad = sigma*2;
        padin = [repmat(in(1,:),[pad 1]); in; repmat(in(end,:),[pad 1])];
        ns = size(padin,1);
        gw = gausswin(floor(ns/2)*2+1,sigma);
        gw = gw(1:ns);
        padout = ifft(bsxfun(@times,fft(padin),ifftshift(gw)));
        out = padout(pad+1:end-pad,:);
end
