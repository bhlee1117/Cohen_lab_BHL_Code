nFiles=length(tiffnames)
countup=0;
for iFile = 1:nFiles
    a=imfinfo(tiffnames{iFile},'tif');
    nFrames=length(a);
    if (nFrames>100)
        countup=countup+1;
        iFile,nFrames
        tiffnames(iFile)
    end
end

countup
    