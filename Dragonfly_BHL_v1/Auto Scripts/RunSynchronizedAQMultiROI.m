

roi = DefineMultiROI;
%%
dmd_load_reference_img;
dmd_set_sequence(hDMDApp,'Camera ROIs',roi);
%%
roi = DefineMultiROI;
%%
[cam_r,cam_c] = ind2sub(size(hDragonflyApp.camApp.cameraMaskImg),find(hDragonflyApp.camApp.cameraMaskImg));
cam_offset = [cam_c(1) cam_r(1)]-1;
cam_active_area = [cam_c(end) cam_r(end)]-[cam_c(1) cam_r(1)]+1;
hDragonflyApp.camApp.CenteredCheckBox.Value = 0;

file_name_prefix  = 'M-YQ0201-7_FOV1_d2_5V_sequential_ramp_c';
for ii=2:length(roi)
    dmd_set_sequence(hDMDApp,'Camera ROIs',roi([1 ii]));
    dmd_send_sequence(hDMDApp);
    
    hDragonflyApp.camApp.ActiveareaDropDown.Value = sprintf('%d x %d',round(roi{ii}(3,1:2)-roi{ii}(1,1:2)));
    hDragonflyApp.camApp.updateActiveArea
    hDragonflyApp.camApp.OffsetDropDown.Value = sprintf('%d, %d',round(roi{ii}(1,1:2)+cam_offset));
    hDragonflyApp.camApp.updateOffset
    hDragonflyApp.appendedfoldernameEditField.Value = ...
        sprintf([file_name_prefix '%g'],ii);
    hDragonflyApp.RunSynchronizedAQ([])
    pause(2)
end

hDragonflyApp.camApp.ActiveareaDropDown.Value = sprintf('%d x %d',cam_active_area);
hDragonflyApp.camApp.updateActiveArea
hDragonflyApp.camApp.OffsetDropDown.Value = sprintf('%d, %d',cam_offset);
hDragonflyApp.camApp.updateOffset