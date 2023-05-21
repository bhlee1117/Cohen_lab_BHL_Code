function pat = dragonfly_special_static_pattern(device,roirows,roicols)
pat = zeros(device.height,device.width);
% normal, random dots for focusing
pat(roicols,roirows) = 1;
pat = pat.*(rand(size(pat))<.005);
return
%%
        datadir = 'R:';

        sdir = dir(fullfile(datadir,'S*'));
        sdir = sdir(cell2mat({sdir.isdir}));
        [~,ix] = sort(arrayfun(@(x)x.name,sdir,'uni',false));
        sdir = sdir(ix(end));

        ramslicepath = fullfile(datadir,sdir.name);

        fdir = dir(fullfile(ramslicepath,'FOV*'));
        fdir = fdir(cell2mat({fdir.isdir}));
        [~,ix] = sort(arrayfun(@(x)x.name,fdir,'uni',false));
        fdir = fdir(ix(end));

        ramfovpath = fullfile(ramslicepath,fdir.name);
        
        dirlist = dir(fullfile(ramslicepath,'*dmd_cal*.mat'));
        [~,ix] = sort(arrayfun(@(x)x.name,dirlist,'uni',false));
%         [~, ix] = sort(cell2mat({dirlist.name}));
%         msgbox(num2str(size(ix)))
        load(fullfile(ramslicepath,dirlist(ix(end)).name));

        dirlist = dir(fullfile(ramfovpath,'*click_spots*.mat'));
        [~,ix] = sort(arrayfun(@(x)x.name,dirlist,'uni',false));
%         [~, ix] = sort(cell2mat({dirlist.name}));
        
        load(fullfile(ramfovpath,dirlist(ix(end)).name));
        stimDotsRowColDMD = [lastClickedX lastClickedY lastClickedX*0+1]*stimCalMatrix;
        stimDotsRowColDMD = round(stimDotsRowColDMD);
        pat = zeros(device.height,device.width);
        for it = 1:size(stimDotsRowColDMD,1)
            pat(stimDotsRowColDMD(it,1),stimDotsRowColDMD(it,2)) = 1;
        end
        pat = imdilate(pat,ones(27,243)); 
        
%         pat = zeros(device.height,device.width);
%         pat(roicols,roirows) = 1;

% for stim clicky calibration
% stimDotsRowColDMD = [
%     599   127
%     617   154
%     197   172
%    -36    245
%    -281   121
%    -458   102
%    -436   222
% ];
% stimDotsRowColDMD = bsxfun(@plus,stimDotsRowColDMD,-[0 170]);
% stimDotsRowColDMD = bsxfun(@rdivide,stimDotsRowColDMD,[-2 .5]);
% stimDotsRowColDMD = bsxfun(@plus,stimDotsRowColDMD,[1024 768]/2);
% stimDotsRowColDMD = round(stimDotsRowColDMD);
% for it = 1:size(stimDotsRowColDMD,1)
%     pat(stimDotsRowColDMD(it,1),stimDotsRowColDMD(it,2)) = 1;
% end

% for stim test
% basepath = fullfile('d:','code','dragonfly runtime');
% dirlist = dir(fullfile(basepath,'*click_spots*.mat'));
% [~, ix] = sort(cell2mat({dirlist.datenum}));
% load(fullfile(basepath,dirlist(ix(end)).name));
% stimDotsRowColDMD = [lastClickedX lastClickedY lastClickedX*0+1]*stimCalMatrix;
% stimDotsRowColDMD = round(stimDotsRowColDMD);
% for it = 1:size(stimDotsRowColDMD,1)
%     pat(stimDotsRowColDMD(it,1),stimDotsRowColDMD(it,2)) = 1;
% end
% pat = imdilate(pat,ones(3));


% for spot pattern test
% load('D:\Code\Hadamard runtime\Analysis\111808_stimulation_masks.mat')
% pat = stimulation_masks(:,:,1)';
