function alignedVolt = transformCamera_O2B(Device_Data, tformReg, AvgVoltageImg, AvgGluImg)
% transformCamera_O2B  Warp voltage image(s) into glutamate camera coordinates.
%
% Usage:
%   alignedVolt = transformCamera_O2B(Device_Data, tformReg, AvgVoltageImg, AvgGluImg)
%
% Inputs:
%   Device_Data    - device metadata cell array
%   tformReg       - affine2d transform (Volt -> Glu)
%   AvgVoltageImg  - H_v x W_v       single voltage image, OR
%                    H_v x W_v x N   stack of N voltage images
%   AvgGluImg      - H_g x W_g       reference glutamate image (2D only, used
%                    H_g x W_g x K     for output sizing; averaged if stacked)
%
% Output:
%   alignedVolt    - H_g x W_g x N   each voltage slice warped into Glu space
%                    (H_g x W_g if input was 2D)

    %-- Collapse GluImg to 2D (only needed for output view sizing)
    if ndims(AvgGluImg) == 3
        G = mean(AvgGluImg, 3);
    else
        G = AvgGluImg;
    end

    %-- Keep voltage stack as-is; remember original dimensionality
    V        = double(AvgVoltageImg);
    nSlices  = size(V, 3);          % 1 for a plain 2D image

    %-- Build spatial referencing objects from camera ROI metadata
    Vmeta = Device_Data{1,3};
    Gmeta = Device_Data{1,4};
    Rvolt = refFromROI(size(V(:,:,1)), double(Vmeta.ROI));
    Rglu  = refFromROI(size(G),        double(Gmeta.ROI));

    %-- Invert transform (tformReg maps Volt->Glu; imwarp needs Glu->Volt)
    tformInv = invert(tformReg);

    %-- Warp each voltage slice independently into Glu coordinates
    % Pre-allocate output in Glu image size
    outSize     = Rglu.ImageSize;           % [H_g  W_g]
    alignedVolt = zeros(outSize(1), outSize(2), nSlices, 'like', V);

    for k = 1:nSlices
        alignedVolt(:,:,k) = imwarp(V(:,:,k), Rvolt, tformInv, 'OutputView', Rglu);
    end

    %-- If input was 2D, return 2D
    if nSlices == 1
        alignedVolt = alignedVolt(:,:,1);
    end
end
