function [blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data)
Rfixed = imref2d(size(Device_Data{1, 6}.refimage.img));
inverseTform = invert(Device_Data{1, 6}.tform);
revertedImage = imwarp(Device_Data{1, 6}.Target, inverseTform,'OutputView',Rfixed);

blueDMDimg=imcrop(revertedImage,double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
bluePatt = bwboundaries(blueDMDimg);