n_rep = 1;
blue_V_max_c = 3;
blue_V_max_w = blue_V_max_c*.4;

blue_V = [linspace(blue_V_max_c*.2,blue_V_max_c,5)' ...
        linspace(max(0.06,blue_V_max_w*.2),blue_V_max_w,5)'];
    
pause_time = 60;

file_name_prefix  = sprintf('M-YQ0201-7_FOV7_d1_c1_wt%gs',pause_time);
file_name_suffix = {
                    '_c_only'    
                    '_w_only'
                    };
                        
waveform_functions = {  
                        'dragonflyMultiPatWaveform_c_only_yq'
                        'dragonflyMultiPatWaveform_w_only_yq'
                        };
for kk = 1:n_rep
for ii=1:size(blue_V)
    
    cnt = 1;
    for jj = randperm(2)
        hDragonflyApp.BlueAmpVEditField.Value = blue_V(ii,jj);
        hDragonflyApp.WaveformFunctionDropDown.Value = waveform_functions{jj};
        hDragonflyApp.appendedfoldernameEditField.Value = ...
            sprintf([file_name_prefix file_name_suffix{jj} '_V_inc%g_rep%g'],ii,kk);
        hDragonflyApp.RunSynchronizedAQ([])
        if ii == size(blue_V,1) && cnt == 2 && kk == n_rep
            continue
        end
        pause(pause_time)
        cnt = cnt+1;
    end
%     hDragonflyApp.WaveformEditField.Value = 'dragonflyMultiPatWaveform_c_yq';
%     hDragonflyApp.appendedfoldernameEditField.Value = ...
%                     sprintf([file_name_prefix '_c_rep%g'],ii);
%     hDragonflyApp.RunSynchronizedAQ([])
%     pause(pause_time)
% %                 
%     hDragonflyApp.WaveformEditField.Value = 'dragonflyMultiPatWaveform_cw_yq';
%     hDragonflyApp.appendedfoldernameEditField.Value = ...
%                     sprintf([file_name_prefix '_cw_rep%g'],ii);
%     hDragonflyApp.RunSynchronizedAQ([])
% 	pause(pause_time)
%     
%     hDragonflyApp.WaveformEditField.Value = 'dragonflyMultiPatWaveform_sw_yq';
%     hDragonflyApp.appendedfoldernameEditField.Value = ...
%                     sprintf([file_name_prefix '_sw_rep%g'],ii);
%     hDragonflyApp.RunSynchronizedAQ([])
%     if ii == n_trial
%         continue
%     end
%     pause(pause_time)
end
end