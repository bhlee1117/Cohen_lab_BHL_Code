function [tiffnames,numframes]= convert_batch_tif()
% function out = convert_batch_tif()
% recognizes .sif file
% reads in .tif per frame from batch conversion in same folder as .sif
% creates single .tif (layered)
% deletes set of .tif per frame generated from batch
% Sabrina Leslie, April 2009

siffiles = dir('*.sif');
siffnames = {siffiles(:).name};
tiffnames = cell(1,length(siffnames));

numframes = zeros(1,length(siffnames)); % stores number of frames per file

for i = 1:length(siffnames)
    for j = 0:9999
        tiffile = dir([siffnames{i}(1:end-4),sprintf('_%04d.tif',j)]);
        if isempty(tiffile)
            break;
        end
    end
    numframes(i) = j;

    %toc;

    %% for a single movie, make one tif file

    if numframes(i)>0

       %img1=imread([siffnames{i}(1:end-4),sprintf('_%04d.tif',0)],'tif');
        %[rows, cols] = size(img1);
        %out = uint16(zeros(numframes(i), rows, cols));
        %out(1,:,:)=img1;

        for j = 1:numframes(i)
            out=imread([siffnames{i}(1:end-4),sprintf('_%04d.tif',j-1)],'tif');
            imwrite(out, [siffnames{i}(1:end-4),'.tif'], 'tif', 'Compression', 'none', 'WriteMode', 'append');
        
            delete([siffnames{i}(1:end-4),sprintf('_%04d.tif',j-1)]);
        
        end; %for
        
        else
        
        numframes(i)=1;

    end %if
    
    tiffnames{i}=[siffnames{i}(1:end-4),'.tif'];
    
end; %for

'done!'