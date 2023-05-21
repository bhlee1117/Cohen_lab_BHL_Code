function alp_patterns = dragonflyL1PatternCalPatterns_yq(liveOffset,dmd,registration_data)
    fig = uifigure;
    selection = uiconfirm(fig,'Please select L1 stimulation pattern sequence','','option',...
                {'Center Surround',...
                'Invert Center Surround',...
                'Ramp Center Surround'},...
                'closefcn',@(h,e) close(fig));
    fig = uifigure;
    wait_input = uiconfirm(fig,'Please select desired figure', '',...
    'option',{'finish'},'closefcn',@(h,e) close(fig));

    hold(gca,'on');
    title('Please clicky Center then Surround')
    npts = 1;
    nroi = 1;
    while(npts > 0)

        [xv, yv] = (getline(gca, 'closed'));
        if size(xv,1) < 3  % exit loop if only a line is drawn
        break
        end

        %draw the bounding polygons and label them
        plot(xv, yv, 'Linewidth', 1);

        roi_points{nroi} = [xv, yv]+liveOffset';
        nroi = nroi + 1;
    end

    L = dmd.device.height;
    W = dmd.device.width;
    pat = false(L,W);
    if nroi>1

        dmd_cam_trans = registration_data(1,:)';
        dmd_rot_dil_mat = registration_data(2:3,:)';
        [x, y] = meshgrid(1:W, 1:L);
        dmd_pixel_pos = cellfun(@(xy) round(dmd_rot_dil_mat\(xy'-repmat(dmd_cam_trans,[1 length(xy)]))),roi_points,'uniformoutput',false);

        in = cellfun(@(pixel_pos) inpolygon(x,y,pixel_pos(1,:),pixel_pos(2,:)),dmd_pixel_pos,'uniformoutput',false);
    end


    switch selection
        case 'Center Surround'
            global dmd_frames_1 dmd_frames_2;
            in{1} = reshape(in{1},[numel(in{1}) 1]);
            in{2} = reshape(in{2},[numel(in{2}) 1]);
            mult1 = ones(1,dmd_frames_1*2); mult1(2:2:end)=0;
            mult2 = ones(1,dmd_frames_2*2); mult2(2:2:end)=0;
            pat_out = cat(3,reshape(in{1}*mult1,L,W,[]),reshape(in{2}*mult2,L,W,[]));
        case 'Invert Center Surround'
            pat = cat(3,pat|(xor(in{1},in{2})),pat|in{2});
        case 'Ramp Center Surround'
            pat = cat(3,pat|in{1},pat,pat|in{2});
    end
    
    figure;moviesc(pat_out);

    alp_patterns = alp_logical_to_btd(permute(pat_out,[2 1 3]));

end