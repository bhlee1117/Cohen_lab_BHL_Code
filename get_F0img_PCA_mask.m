function [ampImg D]=get_F0img_PCA_mask(mov_res,mask,Frame2use)
% get_F0img_PCA_mask - Masked version of get_F0img_PCA.
%
% A 2D binary mask restricts PCA to pixels of interest, reducing
% computation. The output ampImg is full-sized (same as mov_res spatial
% dimensions), with non-mask pixels set to NaN.
%
% Usage:
%   ampImg = get_F0img_PCA_mask(mov_res, mask)
%   ampImg = get_F0img_PCA_mask(mov_res, mask, Frame2use)
%
% Inputs:
%   mov_res    - 3D movie array [Y x X x T]
%   mask       - 2D binary mask [Y x X], logical or 0/1
%   Frame2use  - (optional) frame indices to use; default = all frames
%
% Output:
%   ampImg     - 2D image [Y x X], std of PCA-reconstructed movie within
%                mask; NaN outside mask

if nargin < 3
    Frame2use = [1:size(mov_res,3)];
end

mask = logical(mask);
sz = size(mov_res(:,:,Frame2use));   % [Y X T]
nT = sz(3);

% Extract only masked pixels: [nPix_mask x T]
mov_sub = reshape(mov_res(:,:,Frame2use), [], nT);  % [Y*X x T]
mov_sub = mov_sub(mask(:), :);                       % [nPix_mask x T]

fprintf('Masked PCA: using %d / %d pixels (%.1f%% of frame)\n', ...
    sum(mask(:)), sz(1)*sz(2), 100*sum(mask(:))/(sz(1)*sz(2)));

% PCA on masked pixels only
[u, D, v] = get_eigvector(mov_sub');   % [T x nPix_mask] -> eigenvectors
reshape_u = zeros(sz(1)*sz(2),20);
reshape_u(mask(:),:)= v(:,1:20);
reshape_u = reshape(reshape_u , sz(1), sz(2), []);

figure(999); clf;
imshow2_patch(reshape_u);
n = input("PCs to keep\n");
close(figure(999));

% Print variance explained by the selected PCs
var_total_top20 = sum(D);
var_kept = sum(D(n));
pct_of_top20 = 100 * var_kept / var_total_top20;
fprintf('  Variance explained by selected PCs: %.2f%% \n', pct_of_top20);
fprintf('  Singular values of selected PCs: %s\n', num2str(D(n)'./var_total_top20, '%.4f '));

% PCA reconstruction on masked pixels only

mov_sub_filt=pcafilt(permute(mov_sub,[1 3 2]),n);

% Output: std image, NaN outside mask
ampImg =zeros(sz(1)*sz(2),1);
ampImg(mask(:),1)= std(mov_sub_filt, 0, 3, 'omitnan');
ampImg=reshape(ampImg,sz(1),sz(2),[]);

end
