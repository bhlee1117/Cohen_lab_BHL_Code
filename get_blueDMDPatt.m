function [blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data,mode)
if nargin<2
    mode='normal';
end
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
switch mode
    case 'normal'
        revertedImage = imwarp(double(Device_Data{1, 6}.Target), inverseTform,'OutputView',Rfixed);
    case 'stack'
        Pattern=double(Device_Data{1, 6}.pattern_stack)==1;
        revertedImage = imwarp(Pattern, inverseTform,'OutputView',Rfixed);
end

if sum(sum(double(Device_Data{1, 6}.Target)))>0
    if cam==2
        for z=1:size(revertedImage,3)
            blueDMDimg(:,:,z)=imcrop(revertedImage(:,:,z),double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
        end
    else
        blueDMDimg=imcrop(revertedImage,double(Device_Data{1, 4}.ROI([1 3 2 4]))+[0 0 -1 -1]);
    end
    for z=1:size(revertedImage,3)
        bluePatt(z,:) = bwboundaries(blueDMDimg(:,:,z));
    end

else
    bluePatt=[];     blueDMDimg=[];
end