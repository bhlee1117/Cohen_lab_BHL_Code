n_trial = 1;

pause_time = 60;

file_name_prefix  = 'M-YQ0201-7_FOV1_d1';
file_name_suffix = {
                    '_c'    
                    '_cw'
                    '_sw'
                    };
                        
waveform_functions = {  
                        'dragonflyMultiPatWaveform_c_2color_yq'
                        'dragonflyMultiPatWaveform_cw_2color_yq'
                        'dragonflyMultiPatWaveform_sw_2color_yq'
                        };
                    
dmd_masks = hDMDApp.mask_sequence;
dmd_mask_rois = hDMDApp.mask_sequence_rois;
dmd_mask_c = dmd_masks(:,:,[1 2 2]);
dmd_mask_cw = dmd_masks(:,:,[1 2 3]);
dmd_mask_sw = dmd_masks(:,:,[1 4 3]);
dmd_mask_rois_c = dmd_mask_rois([1 2 2]);
dmd_mask_rois_cw = dmd_mask_rois([1 2 3]);
dmd_mask_rois_sw = dmd_mask_rois([1 4 3]);

dmd_mask_rois_all = {dmd_mask_rois_c,dmd_mask_rois_cw,dmd_mask_rois_sw};
dmd_mask_all = {dmd_mask_c,dmd_mask_cw,dmd_mask_sw};
%%
for ii=1:n_trial   
    
    cnt = 1;
    for jj = randperm(3)
        hDragonflyApp.WaveformFunctionDropDown.Value = waveform_functions{jj};
        hDragonflyApp.appendedfoldernameEditField.Value = ...
            sprintf([file_name_prefix file_name_suffix{jj} '_rep%g'],ii);
        hDMDApp.mask_sequence_rois = dmd_mask_rois_all{jj};
        hDMDApp.mask_sequence = dmd_mask_all{jj};
        dmd_send_sequence(hDMDApp);
        hDragonflyApp.RunSynchronizedAQ([])
        if ii == n_trial && cnt == 3
            continue
        end
        pause(pause_time)
        cnt = cnt+1;
    end

end