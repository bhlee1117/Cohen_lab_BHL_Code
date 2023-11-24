function c_ftprnt=mask_footprint(centers,mov_res,ref_time,rad)
if isempty(ref_time)
    ref_time=[1:size(mov_res,3)];
end

img_cent=zeros(size(mov_res,1),size(mov_res,2),size(centers,1));
for i=1:size(centers,1)
 [x y]=meshgrid([1:size(mov_res,2)],[1:size(mov_res,1)]);
 dist=sqrt(sum(([x(:) y(:)]-centers(i,[1 2])).^2,2));
m=reshape(dist<rad,size(mov_res,1),[]);
mask_cent(i,:)=m(:);
end

D_cent=-bwdist(~(toimg(sum(mask_cent,1),size(mov_res,1),size(mov_res,2))>0));
L=watershed(D_cent); L(~(toimg(sum(mask_cent,1),size(mov_res,1),size(mov_res,2))>0))=0;
nROI = max(L(:));
%imshow2(L,[])


stats2 = regionprops(L, D_cent, 'Area', 'MaxIntensity', 'MinIntensity', 'PixelIdxList');

nROI2 = length(stats2);
for j=1:nROI2; area(j)=stats2(j).Area; end; stats2(find(area<rad^2))=[];
nROI2 = length(stats2);

mov_res_vec = tovec(mov_res(:,:,ref_time));

clear traces c_ftprnt

n_comp = 4; g=1; nFrames=size(mov_res,3); % take 3 PCs per ROI
c_ftprnt = zeros(size(mov_res,1)*size(mov_res,2),nROI2);

for j = 1:nROI2
%     subMov = butterworth_filt(double(mov_r_vec(stats2(j).PixelIdxList,idx_seg)'),4,[10 inf],800)';  % take the pixels in the ROI

    subMov = mov_res_vec(stats2(j).PixelIdxList,:);

    covMat = subMov*subMov';  % PCA within each region
    [V, D] = eig(covMat);
    D = diag(D); 
    D = D(end:-1:1);
    V = V(:,end:-1:1);
    vSign = sign(max(V) - max(-V));  % make the largest value always positive
    V = V.*vSign;

%     eigTraces = V'*subMov;
% figure(8); clf
% for i=1:10
% plot(rescale(eigTraces(i,:)')+i-0.5)
% hold all
% end
% 
% nKeep = 10;
% eigImgs = zeros(size(mov_res,2), size(mov_res,1), nKeep);
% for j = 1:10;
%     eigImgs(:,:,j) = mean(mov_res.*reshape(eigTraces(j,:), [1, 1, size(mov_res,3)]),3);
% end;
% figure; clf;
% for j = 1:nKeep;
%     nexttile([1 1]);
%     imshow2(eigImgs(3:end-3,3:end-3,j), []);
%     title(num2str(j))
% end;

%     mask = mat2gray(sum(V(:,1:n_comp).*D(1:n_comp)'.*(V(:,1:n_comp)>0),2));
    mask = mat2gray(mean(abs(V(:,1:n_comp)).*D(1:n_comp)',2));
    %mask=imgaussfilt(mask,rad/10);
    c_ftprnt(stats2(j).PixelIdxList,j) = mask;
    
%     pause

end
disp(['N = ' num2str(nROI2) ' Cells'])
c_ftprnt = reshape(c_ftprnt,size(mov_res,1),size(mov_res,2),[]);
end