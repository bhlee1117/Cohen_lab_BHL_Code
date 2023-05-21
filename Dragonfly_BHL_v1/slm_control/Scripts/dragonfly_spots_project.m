%%
% S2
% test depth = 1.7969, new depth = 1.7969, nc =  8, p = 144.5
% test depth = 1.7665, new depth = 1.7665, nc = 10, p = 167.3
% test depth = 1.7478, new depth = 1.7478, nc = 10, p =   162
% test depth = 1.7205, new depth = 1.7205, nc =  8, p = 145.2
% 
% load('R:\S2\FOV1\023117_slm_click_spots.mat')
% load('R:\S2\FOV1\023353_slm_click_spots.mat')
% load('R:\S2\FOV1\023550_slm_click_spots.mat')
% load('R:\S2\FOV1\023720_slm_click_spots.mat')

% load('R:\S2\FOV1\024715_optimize_project.mat')
% load('R:\S2\FOV1\031207_optimize_project.mat')
% load('R:\S2\FOV1\031503_optimize_project.mat')
% load('R:\S2\FOV1\032430_optimize_project.mat')

% % 
% % % load 'R:\S1\FOV1\165106_slm_click_spots.mat'
% % % load 'R:\S1\FOV1\165240_slm_click_spots.mat'
% % % load 'R:\S1\FOV1\165407_slm_click_spots.mat'
% % %
% % % load 'R:\S1\FOV1\174213_optimize_project.mat' % D3 z(Arch)=1.73906, 11
% % % load 'R:\S1\FOV1\174816_optimize_project.mat' % D4 z(Arch)=1.72765, 10
% % % load 'R:\S1\FOV1\175151_optimize_project.mat' % D1 z(Arch)=1.75214, 8
% % % 
xyz_offset = [0 0 0];
% xyz_offset = best_xyz_offset;
xyz = slmlocs - xyz_offset;
ttd = [
    0 1 1/3 % 0 sind(60) cosd(60)
    1 0 0
    0 -1/20 1 % cosd(60)/100 sind(60)
    ]*xyz';
% commonz = [
% 0 % Tilt (Y; vertical tilt)
% 0 % Tip (X; horizontal tilt)
% 5 % Oblique astigmatism
% 15*5 % Defocus
% -5*5 % Vertical astigmatism
% 0 % Vertical trefoil
% 0 % Vertical coma
% 0 % Horizontal coma
% 10 % Oblique trefoil
% -10 % Oblique quadrafoil
% 0 % Oblique secondary astigmatism
% 10 % Primary spherical
% -25 % Vertical secondary astigmatism
% -5*10 % Vertical quadrafoil
%     ];

approxSLM = gs(commonz, zpolys, ttd, 10, false, true, false, true, true);

slmph1 = angle(approxSLM);
slmph1 = slmph1*128/pi+128;
slmph1 = max(min(slmph1,255),0);
slmph1 = uint8(slmph1);
% figure windowstyle docked
% imshow(slmph1,[])
% colorbar
projected = true;
prefix = [datestr(now,'HHMMSS') '_slm_click_spots'];
a = dbstack;
if ~isempty(dbstack)
    save(fullfile(ramfovpath,[prefix '.mat']),'slmlocs','slmph1','xyz','ttd','lastClickedX','lastClickedY','priorClickedX','priorClickedY')
end
% calllib('BlinkHdmiSdk', 'Write_image', sdk, slmph1, 1);
fullscreen(rot90(slmph1),3)
