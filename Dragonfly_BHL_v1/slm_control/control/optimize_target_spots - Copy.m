 % VP 2018, YQ2021

function optimize_target_spots(app)

df_app = evalin('base','hDragonflyApp');
daq = evalin('base','hDaq');

%% 
spot_idx = str2num(app.SpotEditField.Value);
if isempty(spot_idx)
    spot_idx = 1:size(app.current_xyz,1);
end

xyz_orig = app.xyz_origin;
current_xyz = app.current_xyz; % spot locations in virtual SLM coordinates
prev_displacement = app.prev_displacement; % previous spot displacement
if isempty(prev_displacement); prev_displacement=ones(size(current_xyz));end

lim_prev_displacement = .1; % if optimized displacement is < 10%, ignore
current_rois = cell2mat(app.current_rois); % spot locations in camera coordinates

a = app.fig_mask_draw.Children;
ref_im = a.Children(end).CData;
[roi_nr,roi_nc] = size(ref_im);

ostep_scale = [1 .4 1]*2; % search space dimensions for [x y z]
osteps = linspace(-1,1,4); % search steps in x y z
oidx = sub2ind([roi_nr roi_nc],round(current_rois(:,2)),round(current_rois(:,1))); % index of clicked camera pixel
nFrames = 10; % N camera frames per step

%%
% step to target more quickly by correcting large registration inaccuracies
% before beginning optimization
switch app.QuickStepCheckBox.Value
    case 1
        delta_lim = 10; % pixels, deviation limit for targeting correction
        for i_spt = spot_idx
            
            cur_frame = wgs_spots(apply_optimal_offsets(current_xyz(i_spt,:)',xyz_orig','apply_trans')',30,1);
            app.SLM.project(cur_frame)
%             pause(.5)
            daq.PololuTrigger = 1;
            pause(.5)
            % take image on current targeting
            temp_frames = zeros(roi_nr,roi_nc,nFrames,'uint16');
            for ii = 1:nFrames
                temp_frames(:,:,ii) = df_app.camApp.snapOneFrame;
            end
            orig_target_data = mean(temp_frames,3);
            orig_target_data([1 end],:) = 0;

            [~,max_idx] = max(orig_target_data(:));
            [max_r,max_c] = ind2sub([roi_nr, roi_nc],max_idx);

            delta_x_cam = current_rois(i_spt,1)-max_c;
            delta_y_cam = current_rois(i_spt,2)-max_r;
            
            if ~any(abs([delta_x_cam delta_y_cam 1])>delta_lim)
                continue
            end
                
            
            mat_rot_dil = app.registration_data(2:4,:)';

            delta_xyz = mat_rot_dil\[delta_x_cam delta_y_cam 1]';
%             idx_correction = abs([delta_x_cam delta_y_cam 1])>delta_lim;
%             current_xyz(i_spt,idx_correction) = current_xyz(i_spt,idx_correction) + delta_xyz(idx_correction);
            current_xyz(i_spt,:) = current_xyz(i_spt,:) + delta_xyz';
            app.current_xyz = current_xyz;
            cur_frame = wgs_spots(apply_optimal_offsets(current_xyz(i_spt,:)',xyz_orig','apply_trans')',30,1);
            app.current_mask = cur_frame;
            app.SLM.project(cur_frame)
%             pause(.5)
            
        end
        cur_frame = wgs_spots(apply_optimal_offsets(current_xyz',xyz_orig','apply_trans')',30,1);
        app.current_mask = cur_frame;
        app.SLM.project(cur_frame)
        
        return
end

%%

% this code block executes one iteration of optimization of SLM pattern to
% maximize static fluorescence from a group of point targets. x,y,z
% location is independently varied and sequentially measured in multiple
% steps. The resulting fluorescence for each cell along these steps is
% independently fitted with a 2nd order polynomial for finding new
% estimates of optimized target in x, y, z. Thus all cells are optimized in
% parallel.
% camera, slm, are initialized earlier.
% this block is run multiple times in a row until convergence is seen in
% intermediate results. 

% alldata = zeros(roi_nr,roi_nc,nFrames,numel(osteps));
avgdata = zeros(roi_nr,roi_nc,numel(ostep_scale),numel(osteps));
allframes = zeros(1152,1920,numel(ostep_scale));
no = numel(osteps);
olist = ceil(((1:no).*(-1).^(1:no))/2) + ceil(no/2); % centered at no/2, right then left at 1 increment
olist = olist(end:-1:1);
for it_sortstep = 1:numel(osteps) % loops over 4-step line search
    it_steps = olist(it_sortstep);
    
    app.SLM.project(zeros(1152,1920))
    for it_coeff = 1:numel(ostep_scale) % loops over x,y,z coordinates separately
        
        % skip the coordinates that are optimized
        % condition empirically determined
        if nnz(abs(prev_displacement(spot_idx,:)) <= lim_prev_displacement) > 1 && ...
                all(abs(prev_displacement(spot_idx,it_coeff)) <= lim_prev_displacement)
            continue
        end
        
        test_xyz_offset = current_xyz;
        
        % add xyz displacement to current xyz
        test_xyz_offset(spot_idx,it_coeff) = current_xyz(spot_idx,it_coeff) + ostep_scale(it_coeff).*osteps(it_steps); 
        allframes(:,:,it_coeff) = wgs_spots(apply_optimal_offsets(test_xyz_offset(spot_idx,:)',xyz_orig','apply_trans')',30,1);
    end
    
    for it_coeff = 1:numel(ostep_scale) % loops over x,y,z coordinates separately
        
        % skip the coordinates that are optimized
        % condition empirically determined
        if nnz(abs(prev_displacement(spot_idx,:)) <= lim_prev_displacement) > 1 && ...
                all(abs(prev_displacement(spot_idx,it_coeff)) <= lim_prev_displacement)
            continue
        end

        app.SLM.project(allframes(:,:,it_coeff))
        pause(.5)

        temp_frames = zeros(roi_nr,roi_nc,nFrames,'uint16');
        for ii = 1:nFrames
            temp_frames(:,:,ii) = df_app.camApp.snapOneFrame;
            pause(.001)
        end
        avgdata(:,:,it_coeff,it_steps) = mean(temp_frames,3);
    end
%%
    fprintf('AQ %d%%\n',ceil(it_sortstep/numel(osteps)*100));

end
%
% fit params
% assume parabolic optim
% find max, optimal scale
% update 
% moviefixsc(squeeze(avgdata(:,:,3,:)))

delta_ostep = zeros(size(current_xyz,1),numel(ostep_scale));
% oidx = oidx(1:2)
for it_coeff = 1:numel(ostep_scale) % loops over x,y,z
    
    % skip the coordinates that are optimized
    % condition empirically determined
    if nnz(abs(prev_displacement(spot_idx,:)) <= lim_prev_displacement) > 1 && ...
                all(abs(prev_displacement(spot_idx,it_coeff)) <= lim_prev_displacement)
            continue
    end
    
    bstepsmov = blur(vm(squeeze(avgdata(:,:,it_coeff,:))),1); % average image for x or y or z
%     bstepsmov(oidx,:)
%     plot(bstepsmov(oidx(1),:)')
    for it_spot = spot_idx
%         plot(osteps,bstepsmov(oidx(it_spot),:))
%         title([it_coeff it_spot])
%         drawnow
%         pause(.1)
        p = polyfit(osteps,bstepsmov(oidx(it_spot),:),2); % parabolic fit of linspace(-1,1,4) vs. measured spot intensity
%         delta_ostep(it_spot,it_coeff) = -p(2)/p(1)/2; % location of max of parabola
        osteps_interp = linspace(-1,1,100);
        [~,loc_max] = max(polyval(p,osteps_interp));
        delta_ostep(it_spot,it_coeff) = osteps_interp(loc_max);%-p(2)/p(1)/2; % location of max of parabola
        
        if range(polyval(p,osteps_interp))/(mean(polyval(p,osteps_interp))-100)<.05
            delta_ostep(it_spot,it_coeff) = 0;
        end
            
    end
    figure(749)
    clf
    subplot(1,3,it_coeff)
    plot(linspace(-1,1,100),polyval(p,linspace(-1,1,100)),osteps,bstepsmov(oidx(it_spot),:))
end
delta_ostep = min(delta_ostep,+1); % limit delta_ostep in: [-1 1]
delta_ostep = max(delta_ostep,-1);
disp ' '
disp(delta_ostep)
current_xyz(spot_idx,:) = current_xyz(spot_idx,:) + delta_ostep(spot_idx,:).*ostep_scale;
prev_displacement(spot_idx,:) = delta_ostep(spot_idx,:);

app.current_xyz = current_xyz;
app.prev_displacement = prev_displacement;
pat = wgs_spots(apply_optimal_offsets(current_xyz',xyz_orig','apply_trans')',30,1);
app.current_mask = pat;
app.SLM.project(pat)
% daq.PololuTrigger = 1;
%%
% figure(748)
% moviefixsc(squeeze(avgdata(:,:,1,:)))
% hold on
% cellfun(@(xyz) plot(xyz(:,1),xyz(:,2),'o'),app.current_rois)
%%
