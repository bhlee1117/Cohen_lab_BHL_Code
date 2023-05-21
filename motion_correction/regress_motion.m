function mov_res = regress_motion(mov,shift_x,shift_y)
% do some simple segmentation    
    gSig = 10; 
    gSiz = 15; 
    psf = fspecial('gaussian', round(gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    Y = imfilter(single(squeeze(mean(mov,3))),psf,'same');
    Y2 = zeros(size(Y));
    Y2(11:end-10,11:end-10) = Y(11:end-10,11:end-10);
    [~,L] = bwboundaries(imbinarize(mat2gray(Y2)));
    pix = regionprops(L, Y2,'PixelIdxList');
    
% do motion regression on segmented regions    
mov_res = tovec(mov);
mov_res_cell = arrayfun(@(x) mov_res(x.PixelIdxList,:),pix,'uniformoutput',false);

shift_x = tovec(shift_x);
shift_y = tovec(shift_y);
shift_x_cell = arrayfun(@(x) mean(shift_x(x.PixelIdxList,:)),pix,'uniformoutput',false);
shift_y_cell = arrayfun(@(x) mean(shift_y(x.PixelIdxList,:)),pix,'uniformoutput',false);

% mov_res_cell = cellfun(@(mov,x,y,xy,xx,yy) SeeResiduals_vec(mov,[x;y],1,0) + mean(mov,2),...
%     mov_res_cell,shift_x_cell, shift_y_cell,...
%     'uniformoutput',false);
reg_coeffs_cell = cellfun(@(mov,x,y) SeeResiduals_vec(movmean(mov-mean(mov,2),50,2),[x;y],0,1),...
    mov_res_cell,shift_x_cell, shift_y_cell,...
    'uniformoutput',false);
mov_res_cell = cellfun(@(mov,x,y,c) mov-c*[x;y],...
    mov_res_cell,shift_x_cell, shift_y_cell, reg_coeffs_cell,...
    'uniformoutput',false);
    for ii = 1:length(pix)
        mov_res(pix(ii).PixelIdxList,:) = mov_res_cell{ii};
    end
    mov_res = toimg(mov_res,size(mov,1),size(mov,2));
end