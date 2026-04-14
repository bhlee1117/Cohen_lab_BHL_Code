function bin2tiff(fpath)

load(fullfile(fpath,"output_data.mat"))
try
    sz=double(Device_Data{1, 2}.ROI([2 4]));
catch
sz=double(Device_Data{1, 3}.ROI([2 4]));
end

mov_mc=double(readBinMov([fpath '/frames1.bin'],sz(2),sz(1)));
write_tif_stack(mov_mc, [fpath '/frames1.tiff'])

end
