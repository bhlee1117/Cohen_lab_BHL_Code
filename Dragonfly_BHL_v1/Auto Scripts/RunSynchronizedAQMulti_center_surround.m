n_trial = 1;

pause_time = 60;

file_name_prefix  = 'M-YQ0201-6_FOV6_d1_c4_.25w_.4s_3V_t150ms_wt60s';
file_name_suffix = {
                    '_c'    
                    '_cw'
                    '_sw'
                    };
                        
waveform_functions = {  
                        'dragonflyMultiPatWaveform_c_yq'
                        'dragonflyMultiPatWaveform_cw_yq'
                        'dragonflyMultiPatWaveform_sw_yq'
                        };
 
for ii=1:n_trial   
    
    cnt = 1;
    for jj = randperm(3)
        hDragonflyApp.WaveformFunctionDropDown.Value = waveform_functions{jj};
        hDragonflyApp.appendedfoldernameEditField.Value = ...
            sprintf([file_name_prefix file_name_suffix{jj} '_rep%g'],ii);
        hDragonflyApp.RunSynchronizedAQ([])
        if ii == n_trial && cnt == 3
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