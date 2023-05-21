%% display img and clicked
dirlist = dir(fullfile(ramfovpath,'*.tif'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
plimg = imread(fullfile(ramfovpath,dirlist(ix(end)).name));
[plr, plc] = size(plimg);
figure windowstyle docked
imshow(plimg,[])
hold on
plot(lastClickedX-cplc/2+plc/2, lastClickedY-cplr/2+plr/2,'o')
prefix = [datestr(now,'HHMMSS') '_spot'];
saveas(gcf,fullfile(ramfovpath,[prefix '.fig']))
