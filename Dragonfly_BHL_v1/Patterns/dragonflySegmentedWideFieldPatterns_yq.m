function alp_patterns = dragonflySegmentedWideFieldPatterns_yq(liveOffset,liveArea,cell_rois,dmd,registration_data)

        
    L = dmd.device.height;
    W = dmd.device.width;
    pat_out=[];
    
    global seg_wide_field_params
    if isempty(seg_wide_field_params)
        opts.Interpreter = 'tex';
        cell_input = inputdlg('\fontsize{10} Please specify [seg_size active_area]','',1,{''},opts);
        seg_wide_field_params = str2num(cell2mat(cell_input));
    end
    seg_size = seg_wide_field_params(1);
    FOV_size = seg_wide_field_params(2:3);
    id11 = [1:seg_size:FOV_size(1)];
    id12 = id11+seg_size-1; id12(id12>FOV_size(1))=FOV_size(1);
    id12 = unique(id12);id11 = id11(1:length(id12));
    id21 = [1:seg_size:FOV_size(2)];
    id22 = id21+seg_size-1; id22(id22>FOV_size(2))=FOV_size(2);
    id22 = unique(id22);id21 = id21(1:length(id22));
    
    fig = uifigure;
    wait_input = uiconfirm(fig,'Please select desired figure', '',...
    'option',{'finish'},'closefcn',@(h,e) close(fig));

    hold(gca,'on');
    
    nroi=1;
    id1_seq = [1:2:length(id11) 2:2:length(id11)];
    for i=1:length(id11)
        pat = false(L,W);
        id1_seq = mod(id1_seq+1,length(id11));
        id1_seq(id1_seq==0)=length(id11);
        for j=1:length(id21)
            
            xv = [id11(id1_seq(j));id11(id1_seq(j));id12(id1_seq(j));id12(id1_seq(j));id11(id1_seq(j))];
            yv = [id21(j);id22(j);id22(j);id21(j);id21(j)];
            
            %draw the bounding polygons and label them
            plot(xv, yv, 'Linewidth', 1);

            roi_points = [xv, yv]+liveOffset';
%             nroi = nroi + 1;
            
            dmd_cam_trans = registration_data(1,:)';
            dmd_rot_dil_mat = registration_data(2:3,:)';
            [x, y] = meshgrid(1:W, 1:L);
            dmd_pixel_pos = round(dmd_rot_dil_mat\(roi_points'-repmat(dmd_cam_trans,[1 length(roi_points)])));

            in = inpolygon(x,y,dmd_pixel_pos(1,:),dmd_pixel_pos(2,:));
            se = strel('square',round(seg_size*.1));
            pat = pat|imdilate(in,se);
        end
        pat_out = cat(3,pat_out,pat);
    end 
    
    figure;moviesc(pat_out)
    
    alp_patterns = alp_logical_to_btd(permute(pat_out,[2 1 3]));

end

% function 