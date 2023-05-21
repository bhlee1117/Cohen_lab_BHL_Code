clear functions
lab_init_device
addpath 'D:\Code\Dragonfly runtime\Scripts'
if(~exist('sdk','var') || isempty(sdk))
    dragonfly_slm_init
end

clear functions
[roirows, roicols] = dragonfly_reload_roi('activity');

ramslicepath = 'R:\S1';
ramfovpath = 'R:\S1\FOV1';

dragonfly_slm_init
dragonfly_cal_slm_display
dragonfly_spots_add_project