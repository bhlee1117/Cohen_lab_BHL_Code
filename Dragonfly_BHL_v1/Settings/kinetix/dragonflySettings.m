% NI DAQ channel definition
outChannelSettings = {
%     'cameraTrigger',    {'Dev1','Port0/Line31', 'OutputOnly'}
    'cameraTrigger',    {'Dev1','Port0/Line26', 'OutputOnly'}
    'grashopperTrigger',{'Dev1','Port0/Line0', 'OutputOnly'}
    'dmdTrigger',       {'Dev1','Port0/Line8',  'OutputOnly'}
    'blueIntensity',    {'Dev1','ao2',          'Voltage'   }
    'orangeIntensity',  {'Dev1','ao0',          'Voltage'   }
    'redShutter',       {'Dev1','Port0/Line22', 'OutputOnly'}
    'limeShutter',      {'Dev1','Port0/Line23', 'OutputOnly'}
    'blueShutter',      {'Dev1','Port0/Line28', 'OutputOnly'}
    'piezoTrigger',     {'Dev1','ao1',          'voltage'   }
    'IRVoltage',        {'Dev1','ao3',          'voltage'   }
%     'PololuTrigger',    {'Dev1','Port0/Line25', 'OutputOnly'}
%     'ODWheelTrigger',   {'Dev1','Port0/Line0',  'OutputOnly'}
    };

inChannelSettings = {
    'encoderA',         {'Dev1','Port0/Line9',  'InputOnly' }
    'encoderB',         {'Dev1','Port0/Line10', 'InputOnly' }   
};

samplingRate = 1e5; 
% samplingRate = 5e5;
% clockSettings = {'External','Dev1/PFI0','ScanClock'};
clockSettings = {'External','Dev1/PFI15','ScanClock'};
% 561 laser control
comPort = 'COM4';
% AOTF com port
comPort_aotf = 'COM10';

% stepper motors serial numbers
ZmotorSerialNumber = 80855715; % z focus linear actuator
WmotorSerialNumber = 80855721; % filter wheel stepper motor

% App figure positions
dragonfly_pos = [980,150,940,584];
camera_pos = [960,850,822,188];
zmotor_pos = [-500, 752, 500, 300];
wmotor_pos = [-500, 452, 500, 300];
dmd_pos = [-260 1 260 800];
slm_pos = [-560 1 300 700];