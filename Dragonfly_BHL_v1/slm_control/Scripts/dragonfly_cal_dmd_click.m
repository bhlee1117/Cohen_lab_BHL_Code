%% DMD click calibrate
% ramspath = 'D:\Code\Dragonfly runtime';
dirlist = dir(fullfile(ramslicepath,'*FFb1_DMD_Cal.tif'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
cal_dmd_plimg = imread(fullfile(ramslicepath,dirlist(ix(end)).name));
[cal_dmd_im_r, cal_dmd_im_c] = size(cal_dmd_plimg);
figure windowstyle docked
imshow((double(cal_dmd_plimg)),[])
% zoom(2)
[xdmd, ydmd] = getpts_o;
%%
stimDotsRowColDMD = [
    599   127
    617   154
    197   172
   -36    245
   -281   121
   -458   102
   -436   222
];
stimDotsRowColDMD = bsxfun(@plus,stimDotsRowColDMD,-[0 170]);
stimDotsRowColDMD = bsxfun(@rdivide,stimDotsRowColDMD,[-2 .5]*4);
stimDotsRowColDMD = bsxfun(@plus,stimDotsRowColDMD,[1024 768]/2);
stimDotsRowColDMD = round(stimDotsRowColDMD);
plxy = [stimDotsRowColDMD];
[~, idplxy] = sort(plxy(:,1));
clxy = [xdmd ydmd];
[~, idclxy] = sort(clxy(:,1));
plot(plxy(idplxy,2),clxy(idclxy,2),'o')
stimCalMatrix = [clxy(idclxy,:) clxy(idclxy,1)*0+1]\[plxy(idplxy,:) plxy(idplxy,1)*0+1];

prefix = [datestr(now,'HHMMSS') '_dmd_cal'];
save(fullfile(ramslicepath,[prefix '.mat']),'stimCalMatrix','cal_dmd_im_c','cal_dmd_im_r')

