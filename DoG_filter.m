function [bp, dog, bg] = DoG_filter(img, sigmaSmall, sigmaLarge, varargin)
% bandpass_DoG: spatial bandpass filter via Difference-of-Gaussians
% img         : 2D image
% sigmaSmall  : small blur (pixel noise suppression), e.g. 0.7–1.2 px
% sigmaLarge  : large blur (background scale), e.g. 6–20 px
%
% Outputs:
%   dog : raw DoG response (can be negative)
%   bp  : bandpassed, clipped to >=0 and normalized (good for visualization)
%   bg  : estimated broad background (Gaussian blurred with sigmaLarge)

p = inputParser;
p.addParameter('Normalize', false, @(x)islogical(x) && isscalar(x));
p.addParameter('ClipNegative', true, @(x)islogical(x) && isscalar(x));
p.parse(varargin{:});
doNorm = p.Results.Normalize;
clipNeg = p.Results.ClipNegative;
bp=[];
for z=1:size(img,3)
    if mod(z,round(size(img,3)/5))==0
    fprintf('\rDoG processing... %6.2f%%', z/size(img,3)*100);
    end
    img(:,:,z) = single(img(:,:,z));

    gSmall = imgaussfilt(img(:,:,z), sigmaSmall, 'FilterSize', max(3, 2*ceil(3*sigmaSmall)+1));
    bg     = imgaussfilt(img(:,:,z), sigmaLarge, 'FilterSize', max(3, 2*ceil(3*sigmaLarge)+1));

    dog = gSmall - bg;

    if clipNeg
        dog(dog < 0) = 0;
    end
    bp(:,:,z) = dog;

    if doNorm
        bp(:,:,z) = bp(:,:,z) - min(tovec(bp(:,:,z)));
        mx = max(tovec(bp(:,:,z)));
        if mx > 0, bp(:,:,z) = bp(:,:,z) ./ mx; end
    end
end
 fprintf('\n');
 disp('DoG done.')
end