function Vout = ApplyWeights(imgs, corrimg, weightimg, offsetimg)
% Vout = ApplyWeights(imgs, corrimg, weightimg, offsetimg)
% 
% Apply a correlation image (M by N), weight image (M by N) and offset 
% image (M by N) to a movie (M by N by K), to produce an estimated voltage 
% Vout (K by 1).
%
% The offset image is subtracted from each frame.  The pixel weight matrix
% is given by weightimg./corrimg.  The pixel weight matrix is contracted 
% with each frame to give the output.
%
% This function can take the output matrices from extractV and generate an
% estimate of the membrane potential.  Useful for training extractV on one
% dataset (or piece thereof) and then applying the weighting to another
% dataset (or piece thereof).
%
% Adam E. Cohen 22 Jan 2011



[ysize, xsize, L] = size(imgs);
offsetmat = repmat(offsetimg, [1 1 L]);
imgs = imgs - offsetmat;
clear offsetmat
scalemat = repmat(weightimg./corrimg, [1 1 L]);
Vout = squeeze(mean(mean(imgs.*scalemat)));
    