function mask_sequence = dmd_set_sequence(app,opt,roi_points)

L = app.dmd.device.height;
W = app.dmd.device.width;
pat = false(L,W);

switch opt
    case 'Append'
        mask = dmd_draw(app,app.DMDMaskGenButtonGroup.SelectedObject.Text);
        pat = cat(3,app.mask_sequence,mask);
    case 'Delete'
        opts.Interpreter = 'tex';
        cell_input = inputdlg('\fontsize{10} Please specify the indices to delete (fmt: ind1 ind2 ... or ind1:ind2)','',1,{''},opts);
        del_ind = str2num(cell2mat(cell_input));
        
        pat(:,:,del_ind) = [];
    case 'Reorder'
        opts.Interpreter = 'tex';
        cell_input = inputdlg('\fontsize{10} Please specify the new indices (fmt: ind1 ind2 ... or ind1:ind2)','',1,{''},opts);
        new_ind = str2num(cell2mat(cell_input));
        
        pat = pat(:,:,new_ind);
    case 'Camera ROIs'
        cam_offset = [app.HorizontalOffsetEditField.Value app.VerticalOffsetEditField.Value];
        if ~exist('roi_points','var')
            [fname,~]=uigetfile('*.mat');
            load(fname,'roi_points')
        end
        
        roi_points = cellfun(@(xy) xy+cam_offset,roi_points,'uniformoutput',false);
        
        dmd_cam_trans = app.registration_data(1,:)';
        dmd_rot_dil_mat = app.registration_data(2:3,:)';
        [x, y] = meshgrid(1:W, 1:L);
        dmd_pixel_pos = cellfun(@(xy) round(dmd_rot_dil_mat\(xy'-repmat(dmd_cam_trans,[1 length(xy)]))),roi_points,'uniformoutput',false);

        % draw dmd patterns
        in = cellfun(@(pixel_pos) inpolygon(x,y,pixel_pos(1,:),pixel_pos(2,:)),dmd_pixel_pos,'uniformoutput',false);
        pat_temp = false(L,W,length(in));
        for i=1:length(in), pat_temp(:,:,i) = pat_temp(:,:,i)|in{i}; end
        
        pat = pat_temp;
end

mask_sequence = pat;
alp_patterns = alp_logical_to_btd(permute(pat,[2 1 3]));
app.dmd.load_sequence(alp_patterns);

app.MaskNoTotalEditField.Value = size(pat,3);
% figure(app.fig_mask_seq_preview);clf
% moviesc(pat);