function ampImg=get_F0img_PCA(mov_res,Frame2use)
if nargin<2
    Frame2use=[1:size(mov_res,3)];
end

sz=size(mov_res);
[u,s,v] = svds(tovec(mov_res),20);
reshape_u=reshape(u,sz(1),sz(2),[]);
figure(999); clf;
imshow2_patch(reshape_u);
n=input("PCs to keep\n");
close(figure(999));
mov_res=pcafilt(mov_res,n);

% kernel = ones(3,3,3); kernel(2,2,2) = 0;
% kernel = kernel/sum(kernel(:));
% 
% base = mean(mov_res, 3);
% 
% dF = mov_res - base;
% dFFilt = imfilter(dF, kernel, 'replicate');
% covMov=dF.*dFFilt;
% covImg = mean(covMov(:,:,avgFrame), 3);
%ampImg = abs(covImg).^.5;

ampImg=std(mov_res(:,:,Frame2use),0,3,'omitnan');
%imshow2(mat2gray(ampImg), [0 .2])
end