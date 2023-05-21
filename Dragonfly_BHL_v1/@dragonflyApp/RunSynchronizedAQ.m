function RunSynchronizedAQ(app,event)



    if app.ZstackButton.Value
        timeout = 2; % s
        originalPosition = app.GotoSpinner.Value;
        app.GotoSpinner.Editable = 'off';
        app.NSpinner.Editable = 'off';
        app.SetFirstZEditField.Editable = 'off';
        app.SetLastZEditField.Editable = 'off';
        app.SnapButton.Enable = 'off';
        app.RunSynchronizedAQButton.Enable = 'off';
        drawnow
        for it = 1:numel(app.stackZPositions)
            app.GotoSpinner.Value = it;
            app.motorZ.moveTo(app.stackZPositions(app.GotoSpinner.Value),timeout)
            runOneAQ
        end
        app.GotoSpinner.Value = originalPosition;
        app.motorZ.moveTo(app.stackZPositions(originalPosition))
        app.GotoSpinner.Editable = 'on';
        app.NSpinner.Editable = 'on';
        app.SetFirstZEditField.Editable = 'on';
        app.SetLastZEditField.Editable = 'on';
        app.SnapButton.Enable = 'on';
        app.RunSynchronizedAQButton.Enable = 'on';
        app.finishSyncAQ
    else
%         it_num = app.RepeatEditField.Value;
%         t_wait = app.WaitTimesEditField.Value;
%         for it_i = 1:it_num
            runOneAQ
%             if it_i < it_num
%                 pause(t_wait)
%             end
%         end
    end
    app.ZstackButton.Value = false;


function runOneAQ
%             if ~isequal(app.WaveformLamp.Color,app.lampColor(true))
%             end
        try %#okay no catch
            app.tim.delete;
        end 

        if app.ZstackButton.Value
            textit = sprintf('%03d_',it);
        else
            textit = '';
        end
		
		app.camApp.fileparams.dcimg_dir = fullfile('E:',[datestr(now,'HHMMSS_') textit app.appendedfoldernameEditField.Value]);
		app.camApp.RunSynchronizedAQButton.ButtonPushedFcn(app,event)

		app.generateWaveform;        
		niDuration = app.queueWaveform;
		
		
        % if app.ZstackButton.Value

            % app.ni6343.syncStartForeground
            % app.camApp.ClosefileButton.ButtonPushedFcn(app,[]);
        % else
%                     additionalWaitInSeconds = .5;
%                     totalWaitInSeconds = ceil((niDuration+additionalWaitInSeconds)*1e3)*1e-3;
%                     app.ni6343.syncStartBackground
%             totalWaitInSeconds = 0;
%             try
                app.ni6343.syncStartForeground;
%             catch
%                 app.finishSyncAQ;
%             end
%             cmdstr = {
%                 'hDragonflyApp.finishSyncAQ;'
%                 };
%             app.tim = timer(...
%                 'TimerFcn', sprintf('%s',cmdstr{:}),...
%                 'StartDelay',totalWaitInSeconds,'ExecutionMode','singleShot');
%             app.tim.start % launch timer to close the file automatically
%                     why
            app.finishSyncAQ
%         pause(3)
        app.writeSettings(fullfile(app.camApp.fileparams.dcimg_dir,[app.dfQuickSettingsDropDown.Value '.set']));
%         try 
        app.dragonfly_write_settings_mat(fullfile(app.camApp.fileparams.dcimg_dir))
%         end %#ok no catch

end
end