function out = medfilt2nan(A, kernelSize)
%MEDFILT2NAN 2-D median filter that ignores NaN values.
%
%   out = MEDFILT2NAN(A, kernelSize)
%     Applies a median filter of size kernelSize to A, ignoring NaNs.
%     kernelSize = [m n] window size (default [3 3]).
%
% Example:
%   A = [1 NaN 3; 4 5 NaN; 7 8 9];
%   B = medfilt2nan(A,[3 3]);

if nargin<2
    kernelSize = [3 3];
end

% pad with NaNs so window stays full at edges
padm = floor(kernelSize(1)/2);
padn = floor(kernelSize(2)/2);
Apad = padarray(A,[padm padn],NaN,'both');

out = NaN(size(A));
for i = 1:size(A,1)
    for j = 1:size(A,2)
        window = Apad(i:i+kernelSize(1)-1, j:j+kernelSize(2)-1);
        % take median ignoring NaNs
        vals = window(~isnan(window));
        if ~isempty(vals)
            out(i,j) = median(vals);
        else
            out(i,j) = NaN;
        end
    end
end
end
