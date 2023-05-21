function spot_downsampl_mask(app,prct)

cam_offset = [app.HorizontalOffsetEditField.Value app.VerticalOffsetEditField.Value];
dmd_cam_trans = app.registration_data(1,:)';
dmd_rot_dil_mat = app.registration_data(2:3,:)';

mask_new = app.current_mask;
spot_list = regionprops(mask_new,'pixelIdxList');
spot_list = arrayfun(@(x) x.PixelIdxList,spot_list,'uniformoutput',false);
idx_rm = randperm(length(spot_list));
idx_rm(1:round(length(idx_rm)*prct/100))=[];

spot_list = cell2mat(spot_list(idx_rm));
mask_new(spot_list) = false;

mask_bd = cellfun(@(x) flip(x,2),bwboundaries(mask_new,'holes'),'uniformoutput',0);
if length(mask_bd)>1
    mask_bd = [mask_bd';num2cell(nan(2,length(mask_bd))',2)'];
    mask_bd = {cell2mat(reshape(mask_bd,[],1))};
end


mask_roi_new = cellfun(@(xy) (round(dmd_rot_dil_mat*xy'+dmd_cam_trans)-cam_offset')',...
                            mask_bd,'uniformoutput',false)';

app.current_rois = mask_roi_new;
app.current_mask  = mask_new;
app.dmd.project(mask_new);
