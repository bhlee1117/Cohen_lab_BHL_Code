function out = write_tif_stack(Imgs, path, jog)
% function out = write_tif_stack(Imgs, path)
% saves a 3-D stack of images as a tiff stack.
% Imgs is m x m x N, where m is the dimension of each image and N is the
% number of images.
% AEC 10/26/05

disp(path);

% Imgs = double(Imgs);
n = size(Imgs,3);
% mn = min(min(min(Imgs)));
% mx = max(max(max(Imgs)));
% Imgs = 32767*(Imgs - mn)/(mx - mn);
Imgs = uint16(Imgs);
disp(size(Imgs));

plane_size = size(Imgs,1)*size(Imgs,2)*16;
max_tif_size = 3.1e10;
n_blocks = uint16(floor(n/jog*plane_size/max_tif_size))+1;
imgs_per_block = uint16(n/n_blocks);

'writing...'
if n_blocks == 1
    for k = 1:jog:n
        imwrite(Imgs(:,:,k), path, 'Compression', 'none', 'WriteMode', 'append');
    end
else
   [fp, name, ext] = fileparts(path);
    for block=0:n_blocks-1
        block_path = [fp '/' name '_block' char(num2str(block+1)) ext];
        for k = (block*imgs_per_block)+1:jog:min(n, (block+1)*imgs_per_block)
            imwrite(Imgs(:,:,k), block_path, 'Compression', 'none', 'WriteMode', 'append');
        end
    end
end
