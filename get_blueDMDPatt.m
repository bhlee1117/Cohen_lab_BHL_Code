function [blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data)
try
Rfixed = imref2d(size(Device_Data{1, 6}.refimage.img));
catch
    try
        cam=2;
    sensorsize=Device_Data{1, 3}.virtualSensorSize; %fusion
    catch
        cam=1;
    sensorsize=size(Device_Data{1, 3}.testimage); %flash
    end
Rfixed = imref2d([sensorsize sensorsize]);    
end
inverseTform = invert(Device_Data{1, 6}.tform);
revertedImage = imwarp(double(Device_Data{1, 6}.Target), inverseTform,'OutputView',Rfixed);

if cam==2
blueDMDimg=imcrop(revertedImage,double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
else
blueDMDimg=imcrop(revertedImage,double(Device_Data{1, 4}.ROI([1 3 2 4]))+[0 0 -1 -1]);
end
bluePatt = bwboundaries(blueDMDimg);