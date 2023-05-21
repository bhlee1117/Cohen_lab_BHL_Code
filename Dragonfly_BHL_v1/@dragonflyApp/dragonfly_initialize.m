function dragonfly_initialize(app)

%
%   2018-2019 Vicente Parot, 2021 Yitong Qi
%   Cohen Lab - Harvard university

%   Add path

addpath(genpath('Hardware'))
addpath(genpath('Patterns'))
addpath(genpath('Settings'))
addpath(genpath('Waveforms'))
addpath(genpath('Hadamard'))
addpath(genpath('Utilities'))
addpath(genpath('dmd_control'))
addpath(genpath('slm_control'))
            
            
try
    dragonflySettings
catch
    edit dragonflySettings
    app.delete
    error('error running settings file')
end


try app.ni6343 = evalin('base','hDaq'); end %#okay no catch
if isempty(app.ni6343) || ~isvalid(app.ni6343)
    fprintf('initializing NI daq ... \n')
    app.ni6343 = nidaq.getInstance(outChannelSettings,inChannelSettings,samplingRate,clockSettings);
    fprintf([8 'done\n'])
    assignin('base','hDaq',app.ni6343);
end

try app.limeLaser = evalin('base','hLimeLaser'); end %#okay no catch
if isempty(app.limeLaser) || ~isvalid(app.limeLaser)
    try
        app.limeLaser = mpbLaser(comPort); 
    catch me
        warning(me.message)
    end
    assignin('base','hLimeLaser',app.limeLaser);

    app.limeLaser.turnOff
end

% try app.AOTF = evalin('base','hAOTF'); end %#okay no catch
% if isempty(app.AOTF) || ~isvalid(app.AOTF)
%     try
%         app.AOTF = AOTF_GH(comPort_aotf); 
%     catch me
%         warning(me.message)
%     end
%     assignin('base','hAOTF',app.AOTF);
% end

try %#okay no catch
%     app.hMotorsFig = evalin('base','hMotorsFig');
    app.motorZ = evalin('base','hMotorZ');
    app.motorW = evalin('base','hMotorW');
end
if  isempty(app.zMotorFig) || ~isvalid(app.zMotorFig) || ...
    isempty(app.wMotorFig) || ~isvalid(app.wMotorFig) || ...
    isempty(app.motorZ)     || ~isvalid(app.motorZ)     || ...
    isempty(app.motorW)     || ~isvalid(app.motorW)
    fprintf('initializing motors ... \n');
    app.zMotorFig = APTmotor.prepareFigure(zmotor_pos);
    app.wMotorFig = APTmotor.prepareFigure(wmotor_pos);
    app.motorZ = APTmotor   (ZmotorSerialNumber,[0, 0, 500, 300],app.zMotorFig);
    app.motorW = filterWheel(WmotorSerialNumber,[0,   0, 500, 300],app.wMotorFig);
    assignin('base','hZMotorFig',app.zMotorFig)
    assignin('base','hWMotorFig',app.wMotorFig)
    assignin('base','hMotorZ',app.motorZ)
    assignin('base','hMotorW',app.motorW)
    app.wMotorFig.Visible = 'off';
    drawnow
    fprintf([8 'done\n']);
end

try app.ODWheel = evalin('base','hODWheel'); end %#okay no catch
if isempty(app.ODWheel) || ~isvalid(app.ODWheel)
    fprintf('initializing ODWheel ... \n');
    app.ODWheel = thorlabsFilterWheel(app.ni6343);
    fprintf('initializing ODWheel ... done\n');
    assignin('base','hODWheel',app.ODWheel);
end

try app.camApp = evalin('base','hCameraApp'); end %#okay no catch
if isempty(app.camApp) || ~isvalid(app.camApp)
    fprintf('initializing camera ... \n');
    app.camApp = cameraApp;
    fprintf('initializing camera ... done\n');
    assignin('base','hCameraApp',app.camApp);
end

try app.dmdApp = evalin('base','hDMDApp'); end %#okay no catch
if isempty(app.dmdApp) || ~isvalid(app.dmdApp)
    fprintf('initializing DMD ... \n');
    app.dmdApp = dmd_control;
    fprintf('initializing DMD ... done\n');
end

% try app.slmApp = evalin('base','hSLMApp'); end %#okay no catch
% if isempty(app.slmApp) || ~isvalid(app.slmApp)
%     fprintf('initializing SLM ... \n');
%     app.slmApp = slm_control;
%     fprintf('initializing SLM ... done\n');
% end

assignin('base','hDragonflyApp',app);



% run camera if not already running
%             try cameraMex startup, end % implements singleton instance

% store internal properties that update on chage of controls only
%             app.liveCenteredROI = app.validateROI(app.ROIDropDown);
%             app.liveExposureTimeSeconds = app.ExposureTimeMsEditField.Value/1000;            
%             app.waveformCenteredROIInput = app.validateROI(app.ROIDropDown_waveform);
%             app.waveformExposureTimeSecondsInput = app.ExposureTimeMsEditField_waveform.Value/1000;            
%             app.waveformCameraFramesInput = app.camFramesInEditField.Value;            

% load default gui settings from selected quick setting
value = app.dfQuickSettingsDropDown.Value;
fname = [value '.set'];
app.readSettings(fname,true)
% app.updateSpacing
app.logStartup
% app.FilterDropDownValueChanged
app.WaveformLamp.Color = app.lampColor(false);
app.PatternsLamp.Color = app.lampColor(false);

%             app.cur_dir = pwd;
app.HadamardSequenceDropDown.Items = {'[11 3]','[19 4]','[27 6]','[35 10]','[59 09]','[63 14]'};
% app.dmd_pat = false(app.dmd.device.height,app.dmd.device.width);
app.WaveformFunctionDropDown.Items = get_file_list(fullfile(pwd,'Waveforms'),'.m');

app.ActivateWaveformCameraFramesCheckBoxValueChanged;
app.ActivateWaveformExposureTimeCheckBoxValueChanged;
app.AcquisitionOptionsButtonGroupSelectionChanged;

app.dragonflyAppUIFigure.Position = dragonfly_pos;
app.camApp.cameraAppUIFigure.Position = camera_pos;
app.dmdApp.UIFigure.Position = dmd_pos;
app.slmApp.UIFigure.Position = slm_pos;

