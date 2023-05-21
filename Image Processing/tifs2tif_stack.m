function out = tifs2tif_stack(filename, firstimage, lastimage)
% function out = tifs2tif_stack(filename, firstimage, lastimage)
% converts a series of 16-bit tif images into a single 16-bit tif stack, and 
% DELETES THE ORIGINAL IMAGES
% AEC 10/26/05

prev_file = 1;

try    % check whether the output file already exists
    imfinfo([filename '_stack.tif']);
catch
    prev_file = 0;
end;

if prev_file;
    'Output file already exists'
    return;
end;

'reading...'
inp = read_tifs(filename, firstimage, lastimage);

'writing...'
for k = 1:(lastimage-firstimage+1);
    imwrite(squeeze(inp(k,:,:)), [filename '_stack.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append');
end;

'deleting...'
for j = firstimage:lastimage
    if j < 10;
        delete([filename '000' num2str(j) '.tif']);        
    elseif j < 100;
        delete([filename '00' num2str(j) '.tif']);            
    elseif j < 1000;
        delete([filename '0' num2str(j) '.tif']);        
    else
        delete([filename num2str(j) '.tif']);        
    end;
end;