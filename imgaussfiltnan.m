function Bout = imgaussfiltnan(A,sigma,varargin)
%IMGAUSSFILTNAN Gaussian filter that ignores NaN values.
%
%   B = IMGAUSSFILTNAN(A, sigma) applies a Gaussian filter with standard
%   deviation sigma to array A while ignoring NaNs. Works like IMGAUSSFILT
%   but treats NaNs as missing data rather than zeros.
%
%   Additional name/value pairs accepted by IMGAUSSFILT can be passed after
%   sigma (e.g. 'FilterSize',[5 5]).
%
% Example:
%   A = peaks(100); A(30:40,30:40)=NaN;
%   B = imgaussfiltnan(A,2);
%
% Byung-Hun / Pawgoomon, 2025

% create Gaussian kernel using MATLAB's imgaussfilt on an impulse
if isempty(varargin)
    % default filter size based on sigma
    ksize = 2*ceil(2*sigma)+1;
else
    % let imgaussfilt handle the parsing but we need kernel explicitly
    % easier: use fspecial
    p = inputParser;
    addOptional(p,'FilterSize',2*ceil(2*sigma)+1);
    parse(p,varargin{:});
    ksize = p.Results.FilterSize;
end
G = fspecial('gaussian', ksize, sigma);

% mask of valid data
Bout=NaN(size(A));
for z=1:size(A,3)
M = ~isnan(A(:,:,z));
Afilled = A(:,:,z);
Afilled(~M) = 0; % replace NaNs with zero for convolution

% convolve both data and mask
num = imfilter(Afilled, G, 'same','replicate');
den = imfilter(double(M), G, 'same','replicate');

B = num ./ den;
B(den==0) = NaN; % places where everything was NaN remain NaN

Bout(:,:,z)=B;
end
end
