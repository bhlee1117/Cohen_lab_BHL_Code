addpath('Auto Scripts')

xy_spinview_taskbar = [1000, 1170];
xy_txtarea = [1563, 87];
xy_start_rec = [1850, 730];
xy_stop_rec = [1750, 730];
% 
% n_trial = 1;
% downsampl = [10:10:50];
% pause_time = 20;




%%
% file_name  = 'M-YQ0201-12_FOV1_20um_mask6_4p7V_6min';
% hDragonflyApp.appendedfoldernameEditField.Value = file_name;
file_name = hDragonflyApp.appendedfoldernameEditField.Value;
mouse_task_bar_click(xy_spinview_taskbar)
paste_to_textarea(xy_txtarea,['E:\' file_name])
mouse_move_click(xy_start_rec)

hDragonflyApp.RunSynchronizedAQ([])
mouse_move_click(xy_stop_rec)