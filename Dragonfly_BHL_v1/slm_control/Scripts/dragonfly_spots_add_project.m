%% select clicky spots
dirlist = dir(fullfile(ramslicepath,'*_dmd_cal.mat'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
load(fullfile(ramslicepath,dirlist(ix(end)).name));

dirlist = dir(fullfile(ramslicepath,'*_slm_cal.mat'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
load(fullfile(ramslicepath,dirlist(ix(end)).name));

dirlist = dir(fullfile(ramfovpath,'*GFP.tif'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
click_plimg = imread(fullfile(ramfovpath,dirlist(ix(end)).name));
[cplr, cplc] = size(click_plimg);
% clear('lastClickedX','lastClickedY')
if exist('lastClickedX','var')
    if projected
        priorClickedX = [priorClickedX; lastClickedX];
        priorClickedY = [priorClickedY; lastClickedY];
    end
else
    priorClickedX = [];
    priorClickedY = [];
end
projected = false;
figure windowstyle docked
imagesc(click_plimg)
colormap gray
axis image
% zoom(2)
hold on
plot(cplc/2+512*[-1 -1 1 1 -1],cplr/2+72*[-1 1 1 -1 -1])
plot(priorClickedX,priorClickedY,'o')
[lastClickedX, lastClickedY] = getpts_o;
slmlocs = [-(lastClickedX+cal_slm_im_c/2-cplc/2) -(lastClickedY+cal_slm_im_r/2-cplr/2) lastClickedX*0+1]*excCalMatrix;
plot(lastClickedX,lastClickedY,'o')
prefix = [datestr(now,'HHMMSS') '_spot_clicky'];
saveas(gcf,fullfile(ramfovpath,[prefix '.fig']))

%% project spots, can be executed by itself also
dragonfly_spots_project