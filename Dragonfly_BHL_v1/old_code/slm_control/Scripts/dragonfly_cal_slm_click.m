%% SLM click calibrate
% basepath = 'D:\Code\Dragonfly runtime';
% ramspath = 'D:\Code\Dragonfly runtime';
dirlist = dir(fullfile(ramslicepath,'*FFb1_SLM_Cal.tif'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
cal_slm_plimg = imread(fullfile(ramslicepath,dirlist(ix(end)).name));
[cal_slm_im_r, cal_slm_im_c] = size(cal_slm_plimg);
figure windowstyle docked
imshow((double(cal_slm_plimg)),[])

% zoom(2)
[xslm, yslm] = getpts_o;
%%
plxy = [excDotsRowColSLM(1:2,:)'];
[~, idplxy] = sort(plxy(:,1));
clxy = [-xslm -yslm];
[~, idclxy] = sort(clxy(:,1));
plot(plxy(idplxy,2),clxy(idclxy,2),'o')
excCalMatrix = [clxy(idclxy,:) clxy(idclxy,1)*0+1]\[plxy(idplxy,:) plxy(idplxy,1)*0+1];

prefix = [datestr(now,'HHMMSS') '_slm_cal'];
save(fullfile(ramslicepath,[prefix '.mat']),'excCalMatrix','commonz','cal_slm_im_r','cal_slm_im_c')
