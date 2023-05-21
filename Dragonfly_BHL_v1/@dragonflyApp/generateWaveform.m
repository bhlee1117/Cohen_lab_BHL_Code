function allWaves = generateWaveform(app)


fprintf('generating waveforms ... \n')

if app.ActivateWaveformExposureTimeCheckBox.Value
    exposureTime = app.waveformExposureTimeSecondsInput;
else
    exposureTime = app.liveExposureTimeSeconds;
end
if app.ActivateWaveformROICheckBox.Value
    centeredROI = app.waveformCenteredROIInput;
else
    centeredROI = app.liveCenteredROI;
end

% nCells = app.nCellsEditField.Value;
% maxBlueIntensity = app.BlueAmpVEditField.Value;
% clear(app.WaveformEditField.Value)
% waveFun = str2func(app.WaveformEditField.Value);
waveFun = str2func(app.WaveformFunctionDropDown.Value);
% recordFrames = app.camFramesInEditField.Value;
allWaves = waveFun(app);
% selectedButton = app.AcquisitionOptionsButtonGroup.SelectedObject;
% text = selectedButton.Text;
% switch text
%     case 'Regular'
%         recordFrames = app.camFramesInEditField.Value;
% %                     allWaves = waveFun(nCells,recordFrames,icolor,maxBlueIntensity,centeredROI,app.camApp.ReadouttimemsEditField.Value);
%         allWaves = waveFun(app);
%     case 'Hadamard'
%         m_q = str2num(app.HadamardSequenceDropDown.Value);
%         recordFrames = m_q(1)+1;
%         sline = 1024;
%         app.FrametimemsEditField.Value = exposureTime+sline*10e-3;
%         allWaves = waveFun(nCells,recordFrames,icolor,maxBlueIntensity,centeredROI,exposureTime);
%     case 'L1 Pattern'
%         recordFrames = app.camFramesInEditField.Value;
%         allWaves = waveFun(recordFrames,icolor,maxBlueIntensity,exposureTime);
% end

app.wavesToSave = allWaves;
app.wavesData = allWaves.amplitude';
app.camPicturesOutEditField.Value = nnz(diff([0 allWaves.amplitude(1,:)])>0);
app.dmdTrigsOutEditField.Value = nnz(diff([0 allWaves.amplitude(2,:)])>0);
try
app.camApp.RecordframesEditField.Value = app.camPicturesOutEditField.Value;
catch
    disp('Camera closed')
end
app.WaveformLamp.Color = app.lampColor(true);     
fprintf('done\n')