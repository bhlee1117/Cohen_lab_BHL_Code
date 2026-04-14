function [alignedGlu, Rglu, Rvolt] = transformCamera(Device_Data, tformReg, AvgVoltageImg, AvgGluImg)
%ALIGNGLUTOVOLTAGE_NOPAD Warp glutamate into voltage coordinates without
% allocating virtual-sensor padded images.
%
% Output is in the VOLTAGE ROI coordinate frame (same size as AvgVoltageImg).

    % Average stacks if needed
    if ndims(AvgVoltageImg) == 3, V = mean(AvgVoltageImg, 3); else, V = AvgVoltageImg; end
    if ndims(AvgGluImg)     == 3, G = mean(AvgGluImg, 3);     else, G = AvgGluImg;     end

    Vmeta = Device_Data{1,3}; % voltage info: ROI = [x y w h]
    Gmeta = Device_Data{1,4}; % glu info

    % Build spatial referencing for CROPPED images placed in "virtual sensor" world coords
    Rglu  = refFromROI(size(G), double(Gmeta.ROI));
    Rvolt = refFromROI(size(V), double(Vmeta.ROI));

    % Warp glu into voltage coordinates (but only over the voltage ROI area)
    alignedGlu = imwarp(G, Rglu, tformReg, 'OutputView', Rvolt);
end


function R = refFromROI(imgSize, ROI)
% ROI assumed [x y w h] with 0-based-like offsets as in your padding code.

    x = ROI(1);  y = ROI(3);
    w = ROI(2);  h = ROI(4);

    if imgSize(2) ~= w || imgSize(1) ~= h
        error('Image size (%dx%d) does not match ROI (%dx%d).', imgSize(1), imgSize(2), h, w);
    end

    xWorld = [x + 0.5, x + w + 0.5];
    yWorld = [y + 0.5, y + h + 0.5];

    R = imref2d([h w], xWorld, yWorld);
end
