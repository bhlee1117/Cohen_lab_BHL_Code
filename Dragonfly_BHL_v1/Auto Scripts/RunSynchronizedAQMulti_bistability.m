addpath('Auto Scripts')

xy_spinview_taskbar = [1000, 1170];
xy_txtarea = [1548, 80];
xy_start_rec = [1850, 730];
xy_stop_rec = [1750, 730];

n_trial = 5;
downsampl = [10:20:50];
pause_time = 40;


file_name_prefix = 'M-YQ-0201-12_FOV1_20um_1V_downsampl';

%%
for i_trial = 3:n_trial
    for ii = downsampl

        for jj = 1:round(1/(ii*1e-2))
            file_name = sprintf([file_name_prefix num2str(ii) '_rpt' num2str(jj) '_trial' num2str(i_trial)]);
            hDragonflyApp.appendedfoldernameEditField.Value = file_name;

            hDMDApp.FunArgsTextArea.Value = [num2str(ii) ',[1 2],' num2str(jj)];
            dmd_run_quick_function(hDMDApp);
            dmd_send_sequence(hDMDApp);

            mouse_task_bar_click(xy_spinview_taskbar)
            paste_to_textarea(xy_txtarea,['E:\' file_name])
            mouse_move_click(xy_start_rec)

            hDragonflyApp.RunSynchronizedAQ([])
            mouse_move_click(xy_stop_rec)
            figure(90);clf
            moviefixsc(hDMDApp.mask_sequence)
            pause(pause_time)
        end

    end
end