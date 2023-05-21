function out = toimg(mat, nr, nc);

% function out = toimg(mat, nr, nc);
% converts a matrix of size [nr*nc, nt] to a movie of size [nr, nc, nt];
%
% modified from code by Vicente Parot
if ~exist('nc','var') && numel(nr)==2
    nc = nr(2);
    nr = nr(1);
end
    
out = reshape(mat, nr, nc, []);