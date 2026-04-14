
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