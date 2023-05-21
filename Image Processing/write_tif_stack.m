function out = write_tif_stack(Imgs, path)
% function out = write_tif_stack(Imgs, path)
% saves a 3-D stack of images as a tiff stack.
% Imgs is m x m x N, where m is the dimension of each image and N is the
% number of images.
% AEC 10/26/05


n = size(Imgs,3);

mn = min(min(min(Imgs)));
mx = max(max(max(Imgs)));
Imgs = 32767*(Imgs - mn)/(mx - mn);
Imgs = uint16(Imgs);

'writing...'
for k = 1:n;
    imwrite(Imgs(:,:,k), path, 'Compression', 'none', 'WriteMode', 'append');
end;
