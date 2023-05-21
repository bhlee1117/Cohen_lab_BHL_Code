function dmd_arb_operation(app,idx)

se = strel('disk',4,8);
cam_offset = [app.HorizontalOffsetEditField.Value app.VerticalOffsetEditField.Value];
dmd_cam_trans = app.registration_data(1,:)';
dmd_rot_dil_mat = app.registration_data(2:3,:)';
% mask_sequence_new = xor(imdilate(app.mask_sequence(:,:,idx(1)),se),app.mask_sequence(:,:,idx(2)));
mask_sequence_new = and(app.mask_sequence(:,:,idx(1)),app.mask_sequence(:,:,idx(2)));

mask_bd = cellfun(@(x) flip(x,2),bwboundaries(mask_sequence_new,'holes'),'uniformoutput',0);
if length(mask_bd)>1
    mask_bd = [mask_bd';num2cell(nan(2,length(mask_bd))',2)'];
    mask_bd = {cell2mat(reshape(mask_bd,[],1))};
end


mask_roi_new = cellfun(@(xy) (round(dmd_rot_dil_mat*xy'+dmd_cam_trans)-cam_offset')',...
                            mask_bd,'uniformoutput',false)';

app.mask_sequence_rois(end+1) = mask_roi_new;
app.mask_sequence(:,:,end+1)  = mask_sequence_new;
app.MaskNoTotalEditField.Value = size(app.mask_sequence,3);