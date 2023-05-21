function out = write_jpgs(data, filename, firstimage, lastimage)
% function out = write_jpgs(data, filename, firstimage, lastimage)
% writes images from a 3-d array into a series of jpeg files
% AEC 4/19/05



for j = firstimage:lastimage
    if j < 10;
        imwrite(squeeze(data(j,:,:)), [filename '000' num2str(j) '.tif'], 'Compression', 'none');        
    elseif j < 100;
        imwrite(squeeze(data(j,:,:)), [filename '00' num2str(j) '.tif'], 'Compression', 'none');     
    elseif j < 1000;
        imwrite(squeeze(data(j,:,:)), [filename '0' num2str(j) '.tif'], 'Compression', 'none');     
    else
        imwrite(squeeze(data(j,:,:)), [filename num2str(j) '.tif'], 'Compression', 'none');     
    end;
end;
