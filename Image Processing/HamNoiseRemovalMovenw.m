% Hamamatsu camera noise removal wrapper. Handles noise removal from raw
% movies with multiple recording sections
% Eli N. Weinstein, Cohen Lab, 5/5/2015

% Input
% ims -> Input movie (should have the 100 or 400 baseline counts
% subtracted, but no other corrections). Number of pixels by number of
% frames
% varargin{1} (secs) -> cell array consisting of 1 by 2 vectors describing the
% starting and ending frames of each section (continuous camera
% recording) within the larger trace
% varargin{2} (satcoeff) -> Saturation of noise frames (exposure of noise
% frames/exposure of non-noise frames)
% Output
% imsc -> Corrected movie
function imsc = HamNoiseRemovalMovenw(ims, varargin)
%%
if numel(varargin) > 0 && ~isempty(varargin{1})
    secs = varargin{1};
else
    secs = {[1, size(ims, 2)]};
end
if numel(varargin) > 1 && ~isempty(varargin{2})
    osat = 1; % whether to use outside or measured satcoeff
    satcoeff = varargin{2};
else
    osat = 0;
end
%% Filter movie to maximize noise
% Calculate lag-4 autocorrelation; there's a slight error here from
% the different sections, but it should hardly matter
lag4img = sqrt(sum(ims(:, 1:(size(ims, 2) - 5 + 1)).*ims(:, 5:size(ims, 2)), 2)./sum(ims.*ims, 2));
lag4img = lag4img./sum(sum(lag4img)); % normalize

lag4tr  = sum(ims.*repmat(lag4img, 1, size(ims, 2)), 1);
%% Call main noise removal function, filter movie sections
imsc = zeros(size(ims));
for kk = 1:numel(secs)
    thsto = lag4tr(secs{kk}(1):secs{kk}(2));
    [~, noiset, satcoeffm] = HamNoiseRemoval500enw(thsto);
    if osat == 0
        satcoeff = satcoeffm;
    end
    imsc(:, secs{kk}(1):secs{kk}(2)) = ims(:, secs{kk}(1):secs{kk}(2))...
        .*repmat(~noiset + (noiset./satcoeff), size(ims, 1), 1);
end
    