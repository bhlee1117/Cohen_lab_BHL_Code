function ampImg=get_F0img_PCA(mov_res,Frame2use)
if nargin<2
    Frame2use=[1:size(mov_res,3)];
end

sz=size(mov_res(:,:,Frame2use));
%[u,s,v] = svds(tovec(mov_res(:,:,Frame2use)),20);
[u D v]=get_eigvector(tovec(mov_res(:,:,Frame2use))');
reshape_u=reshape(v(:,1:20),sz(1),sz(2),[]);
figure(999); clf;
imshow2_patch(reshape_u);
n=input("PCs to keep\n");
close(figure(999));

% Print variance explained by the selected PCs
var_each = D;                          
var_total_top20 = sum(D);       % total variance captured by top 20 PCs
var_kept = sum(var_each(n));         % variance captured by selected n PCs
pct_of_top20 = 100 * var_kept / var_total_top20;
fprintf('  Variance explained by selected PCs: %.2f%% (of top-20 variance)\n', pct_of_top20);
fprintf('  Singular values of selected PCs: %s\n', num2str(D(n)'./var_total_top20, '%.4f '));

mov_res=pcafilt(mov_res(:,:,Frame2use),n);

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

ampImg=std(mov_res,0,3,'omitnan');
%imshow2(mat2gray(ampImg), [0 .2])
end