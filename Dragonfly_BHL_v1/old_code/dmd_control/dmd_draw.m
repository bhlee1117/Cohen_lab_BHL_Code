function mask = dmd_draw(app,opt)

if isvalid(app.fig_mask_draw)
    figure(app.fig_mask_draw)
else
    app.fig_mask_draw = load_reference_img();
end

L = app.dmd.device.height;
W = app.dmd.device.width;
pat = false(L,W);

cam_offset = [app.HorizontalOffsetEditField.Value app.VerticalOffsetEditField.Value];

switch opt
    case 'Arb Shape'
        
        % start cliky loop
        hold on
        title('Please clicky ROIs')
        npts = 1;
        nroi = 1;
        while(npts > 0)

            [xv, yv] = (getline(gca, 'closed'));
            if size(xv,1) < 3  % exit loop if only a line is drawn
            break
            end

            %draw the bounding polygons and label them
            plot(xv, yv, 'Linewidth', 1);

            roi_points{nroi} = [xv, yv]+cam_offset;
            nroi = nroi + 1;
        end

        
        
    case 'Rectangle'
        
        try %#okay uncaught
            myrect = getrect;   
        end
        if exist('myrect','var') && isnumeric(myrect) && numel(myrect)==4
            roi_points = {myrect(1:2)+myrect(3:4).*[0 0 1 1 0; 0 1 1 0 0]'+cam_offset};
            hold on
            plot(roi_points{1}(:,1),roi_points{1}(:,2))
            hold off
        end
        
    case 'Circle'
    case 'Load Picture'
end


try %#okay uncaught
% convert camera coordinate to dmd coordinate
dmd_cam_trans = app.registration_data(1,:)';
dmd_rot_dil_mat = app.registration_data(2:3,:)';
[x, y] = meshgrid(1:W, 1:L);
dmd_pixel_pos = cellfun(@(xy) round(dmd_rot_dil_mat\(xy'-repmat(dmd_cam_trans,[1 length(xy)]))),roi_points,'uniformoutput',false);

% draw dmd patterns
in = cellfun(@(pixel_pos) inpolygon(x,y,pixel_pos(1,:),pixel_pos(2,:)),dmd_pixel_pos,'uniformoutput',false);
for i=1:length(in), pat = pat|in{i}; end

app.dmd.project(double(pat))
mask = pat;
catch
    disp('Pattern calculation error')
    mask = 0;
end