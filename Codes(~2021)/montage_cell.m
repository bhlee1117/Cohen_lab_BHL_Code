clear
fileFolder = fullfile('E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Classification\OE6\Ret\Non_TXN');
imds = imageDatastore(fileFolder,'IncludeSubfolders',true,'LabelSource',...
    'foldernames');

montage(imds, 'Size', [40 10],'DisplayRange',[100 1000]);