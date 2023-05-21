function spot_downsampl_sequence(app,prct,idx,rpt)
cam_offset = [app.HorizontalOffsetEditField.Value app.VerticalOffsetEditField.Value];
dmd_cam_trans = app.registration_data(1,:)';
dmd_rot_dil_mat = app.registration_data(2:3,:)';

spot_rois = app.mask_sequence_rois{idx(1)};
spot_rois_idx = find(ismissing(spot_rois(:,1)));

mask_new = logical(app.mask_sequence(:,:,idx(1)));
spot_list = regionprops(mask_new,'pixelIdxList');
spot_list = arrayfun(@(x) x.PixelIdxList,spot_list,'uniformoutput',false);

idx_bin = 0:length(spot_list)/(1/prct*100):length(spot_list);
if idx_bin(end)~=length(spot_list), idx_bin(end+1) = length(spot_list);end
if exist('rpt','var'), rng('default');end
idx = 1:length(spot_list);
idx_rm = randperm(length(spot_list));

if exist('rpt','var')
    idx_rm(idx>idx_bin(rpt) & idx<=idx_bin(rpt+1)) = [];
else
    idx_rm(idx>idx_bin(1) & idx<=idx_bin(2)) = [];
end

spot_list = cell2mat(spot_list(idx_rm));
mask_new(spot_list) = false;

spot_rois_new = {};
for ii = 1:length(spot_rois_idx)
    if ii == 1
        spot_rois_new{ii} = spot_rois(1:spot_rois_idx(1),:);
    else
        spot_rois_new{ii} = spot_rois(spot_rois_idx(ii-1)+1:spot_rois_idx(ii),:);
    end
end

spot_rois_new(idx_rm) = [];
spot_rois_new = cell2mat(spot_rois_new');

app.mask_sequence_rois{idx(2)} = spot_rois_new;
app.mask_sequence(:,:,idx(2))  = mask_new;