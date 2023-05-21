function alp_patterns = dragonflyL1PatternPatterns_yq(liveOffset,liveArea,cell_rois,dmd,registration_data)
    fig = uifigure;
    selection = uiconfirm(fig,'Please select L1 stimulation pattern sequence','','option',...
                {'Center Surround',...
                'Invert Center Surround',...
                'Center+Inverted Center Surround',...
                'Center+Inverted Center Surround RNG'},...
                'closefcn',@(h,e) close(fig));
    
    fig = uifigure;
    pattern_input = uiconfirm(fig,'Pattern: predefined or manual?', '',...
    'option',{'Cell Targeted Manual','Manual','Predefined'},'closefcn',@(h,e) close(fig));        
    
    fig = uifigure;
    wait_input = uiconfirm(fig,'Please select desired figure', '',...
    'option',{'finish'},'closefcn',@(h,e) close(fig));

    hold(gca,'on');
    if strcmp(pattern_input,'Cell Targeted Manual')
        cellfun(@(x) plot(x(:,1),x(:,2),'r'),cell_rois)
    end

    switch pattern_input
        case 'Predefined'
            h=gca; im = get(h,'children');im=im(end);
            title('Predefined Pattern')
            roi_points = cell_rois;
            
            i_roi = ceil(rand()*length(roi_points));
            
            fig = uifigure;
            accept = uiconfirm(fig,sprintf('Current ROI: %g. Do you accept?',i_roi),'','option',...
            {'Yes',...
            'No',...
            'Cancel'},...
            'closefcn',@(h,e) close(fig));
        
            while ~strcmp(accept,'Yes')
                if strcmp(accept,'Cancel')
                    return
                else
                    i_roi = ceil(rand()*length(roi_points));
            
                    fig = uifigure;
                    accept = uiconfirm(fig,sprintf('Current ROI: %g. Do you accept?',i_roi),'','option',...
                    {'Yes',...
                    'No',...
                    'Cancel'},...
                    'closefcn',@(h,e) close(fig));
                end
            end
            
            plot(h,roi_points{i_roi}(:,1),roi_points{i_roi}(:,2),'r','linewidth',1)
            [yb,xb]=size(im.CData);
            [roi_center(1),roi_center(2)] = centroid(polyshape(roi_points{i_roi}));
            roi_poly =    nsidedpoly(30,'center',...
                            roi_center,...
                            'radius',...
                            sqrt(polyshape(roi_points{i_roi}).area*5));
            roi_poly = roi_poly.Vertices;
            roi_points={};
            roi_points{1} = [roi_poly;roi_poly(1,:)];
            roi_points{2} = [[1 1 xb xb 1]' [1 yb yb 1 1]'];
            plot(h,roi_points{1}(:,1),roi_points{1}(:,2),roi_points{2}(:,1),roi_points{2}(:,2))
            roi_points = cellfun(@(x) x+liveOffset',roi_points,'uniformoutput',false);
            nroi=length(roi_points);
        otherwise
            title('Please clicky Center')
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
            
            roi_points{end+1}(:,1) = [0 0 liveArea(1) liveArea(1) 0]'+liveOffset(1);
            roi_points{end}(:,2) = [0 liveArea(2) liveArea(2) 0 0]'+liveOffset(2);
            nroi=nroi+1;
            plot(roi_points{end}(:,1)-liveOffset(1), roi_points{end}(:,2)-liveOffset(2), 'Linewidth', 1);
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
        
        if strcmp(pattern_input,'Cell Targeted Manual')
            dmd_pixel_pos_cell = cellfun(@(xy) round(dmd_rot_dil_mat\((xy+liveOffset')'-repmat(dmd_cam_trans,[1 length(xy)]))),cell_rois,'uniformoutput',false);

            pat_cell = cellfun(@(pixel_pos) inpolygon(x,y,pixel_pos(1,:),pixel_pos(2,:)),dmd_pixel_pos_cell,'uniformoutput',false);
            for i=1:length(pat_cell),pat = pat|pat_cell{i};end
            se = strel('disk',5,8);
            pat = imdilate(pat,se);
        end
    end
    
    if strcmp(pattern_input,'Cell Targeted Manual')
        switch selection
            case 'Center Surround'
                pat = cat(3,pat&in{1},pat&in{2});
            case 'Invert Center Surround'
                pat = cat(3,pat&(xor(in{1},in{2})),pat&in{2});
            case 'Ramp Center Surround'
                pat = cat(3,pat&in{1},pat,pat&in{2});
            case 'Center+Inverted Center Surround'
                pat = cat(3,pat&in{1},pat&in{3},pat&(xor(in{2},in{3})),pat&in{3});
            case 'Center+Inverted Center Surround RNG'
                pat = cat(3,pat&in{1},pat&in{3},pat&(xor(in{2},in{3})),pat&in{3});
        end
    else
        switch selection
            case 'Center Surround'
                pat = cat(3,pat|in{1},pat|in{2});
            case 'Invert Center Surround'
                pat = cat(3,pat|(xor(in{1},in{2})),pat|in{2});
            case 'Ramp Center Surround'
                pat = cat(3,pat|in{1},pat,pat|in{2});
            case 'Center+Inverted Center Surround'
                pat = cat(3,pat|in{1},pat|in{3},pat|in{2},pat|in{3});
%                 pat = cat(3,pat|in{1},pat|in{3},pat|(xor(in{2},in{3})),pat|in{3});
            case 'Center+Inverted Center Surround RNG'
                pat = cat(3,pat|in{1},pat|in{3},pat|(xor(in{2},in{3})),pat|in{3});
        end
    end
    
    global input_repeat_wait;
    if ~isempty(input_repeat_wait)
        if contains(selection,'Center+Inverted')
            pat = repmat(pat,[1 1 input_repeat_wait(1)/2]);
        else
            pat = repmat(pat,[1 1 input_repeat_wait(1)]);
        end
        
        if contains(selection,'RNG')
            n_seq = input_repeat_wait(1);
            seq_randperm_idx = randperm(n_seq);
            seq_randperm_idx = reshape([seq_randperm_idx*2-1;seq_randperm_idx*2],1,[]);
            pat = pat(:,:,seq_randperm_idx);
        end
            
    end
    
%     figure;moviesc(pat)

    alp_patterns = alp_logical_to_btd(permute(pat,[2 1 3]));

end