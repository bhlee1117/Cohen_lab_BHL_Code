function pat = dragonfly_display_dmd_cal(device)
pat = zeros(device.height,device.width);
% normal, random dots for focusing
% pat(roicols,roirows) = 1;
% pat = pat.*(rand(size(pat))<.005);

% for stim clicky calibration
stimDotsRowColDMD = [
    599   127
    617   154
    197   172
   -36    245
   -281   121
   -458   102
   -436   222
];
stimDotsRowColDMD = bsxfun(@plus,stimDotsRowColDMD,-[20 210]);
stimDotsRowColDMD = bsxfun(@rdivide,stimDotsRowColDMD,[-2 .5]*4);
stimDotsRowColDMD = bsxfun(@plus,stimDotsRowColDMD,[1024 768]/2);
stimDotsRowColDMD = round(stimDotsRowColDMD);
for it = 1:size(stimDotsRowColDMD,1)
    pat(stimDotsRowColDMD(it,1),stimDotsRowColDMD(it,2)) = 1;
end

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
