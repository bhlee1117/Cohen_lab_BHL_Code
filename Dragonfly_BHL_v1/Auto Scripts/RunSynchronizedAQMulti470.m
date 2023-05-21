
blue_amps = ones(10,1);

if ~strcmp('RFP',hMotorW.currentFilter)
    hMotorW.switchTo('RFP');
%     hDragonflyApp.
end
pause(1)
for ii=1:length(blue_amps)   
    blue_amps_i = blue_amps(ii);
    hDragonflyApp.BlueAmpVEditField.Value = blue_amps_i;
    hDragonflyApp.appendedfoldernameEditField.Value = ...
                    sprintf('YQ0303L_JFx608_DIV16_FOV20_cell13_1V_rep%g',ii);
%                 sprintf('YQ0103_p1_div19_FOV5_%.1gV',blue_amps_i);
                
    hDragonflyApp.RunSynchronizedAQ([])
    if ii == length(blue_amps)
        continue
    end
    pause(20)
end