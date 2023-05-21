function mask = slm_draw(app,opt,type)
%% handle plotting routine

if strcmp(opt,'Load Picture')
    [fname,fpath]=uigetfile('C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\*.*'); if fname==0;return;end
    mask = imread(fullfile(fpath,fname));
    app.SLM.project(mask)
    return
end

if ~isa(app.fig_mask_draw,'double') && isvalid(app.fig_mask_draw)
    figure(app.fig_mask_draw)
else
    app.fig_mask_draw = slm_load_reference_img();
end
Im = findobj(get(app.fig_mask_draw,'children'),'type','Image'); im = Im.CData;

switch opt
    case 'Arb Shape'
        roi_points = click_im(im,'clicky',app.fig_mask_draw);
    case 'Rectangle'
        roi_points = click_im(im,'rect',app.fig_mask_draw);
    case 'Spots'
        roi_points = click_im(im,'pts',app.fig_mask_draw);
        roi_points = mat2cell(roi_points,ones(1,size(roi_points,1)));
    case 'Load Picture'

end

switch type.mode
    case 'Gen'
        app.current_rois = roi_points;
    case 'Add'
        app.current_rois = [app.current_rois; roi_points];
end
hold(app.fig_mask_draw.Children,'off')

cam_offset = [app.HorizontalOffsetEditField.Value app.VerticalOffsetEditField.Value]; % get camera offset
roi_points = cellfun(@(pts) pts+cam_offset,roi_points,'uniformoutput',false); % apply camera offset

%% convert camera coordinate to slm coordinate
xyz_origin = app.xyz_origin;
try 

slm_cam_trans = app.registration_data(1,:)';
slm_rot_dil_mat = app.registration_data(2:4,:)';
slm_pixel_pos = cellfun(@(xy) (slm_rot_dil_mat\([xy 1]'-repmat(slm_cam_trans,[1 size(xy,1)])))',roi_points,'uniformoutput',false);
switch type.type
    case 'Area'
        target = gs_target_gen(slm_pixel_pos,'roi');
        hold(app.fig_mask_draw.Children,'on')
        cellfun(@(pts) plot(app.fig_mask_draw.Children,pts(:,1),pts(:,2)),roi_points)
        pat = gs(target,200);
    case 'Spots'

        switch opt
            case 'Spots'
                xyz = cell2mat(slm_pixel_pos);
                switch type.mode
                    case 'Add'
                        xyz = [app.current_xyz; xyz];
                end
                pat = wgs_spots(apply_optimal_offsets(xyz',xyz_origin','apply_trans')',30,1);
                app.prev_displacement = [];
            otherwise
                xy = gs_target_gen(slm_pixel_pos,'roi2spt');
                z = zeros(size(xy,1),1);
                xyz = [xy z];
                pat = wgs_spots(apply_optimal_offsets(xyz',xyz_origin','apply_trans')',30,1);
                hold(app.fig_mask_draw.Children,'on')
                cellfun(@(pts) plot(app.fig_mask_draw.Children,pts(:,1),pts(:,2),'*'),roi_points)
        end
        app.current_xyz = xyz;
end

% assignin('base','slm_pat',pat)

app.SLM.project(pat)
% slm_project
% calllib('Blink_C_wrapper', 'Write_image', rot90(pat,3), true)
mask = pat;
catch
    disp('Pattern calculation error')
    mask = 0;
end