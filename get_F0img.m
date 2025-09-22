function ampImg=get_F0img(mov_res,avgFrame)
if nargin<2
    avgFrame=[1:size(mov_res,3)];
end

kernel = ones(3,3,3); kernel(2,2,2) = 0;
kernel = kernel/sum(kernel(:));

base = mean(mov_res, 3);

dF = mov_res - base;
dFFilt = imfilter(dF, kernel, 'replicate');
covMov=dF.*dFFilt;
covImg = mean(covMov(:,:,avgFrame), 3);
ampImg = abs(covImg).^.5;
%imshow2(mat2gray(ampImg), [0 .2])
end