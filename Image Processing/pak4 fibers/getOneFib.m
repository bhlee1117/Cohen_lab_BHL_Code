function fiberImg = getOneFib(X0, Y0, dW, img, varargin);
% function fiberImg = getOneFib(X0, Y0, img);
%
% Pull an image of a fiber from an image.
% X0 and Y0 are each 2-element vectors which specify the ends of the fiber.
% The extracted image has 2*dW+1 rows and a boundary of width dW on either
% end of the fiber.
% img can be either black-and-white or color
% varargin: if included, should be an option (OPT) structure with field
% OPT.USE_IMPROFILE (true = use improfile [slower]; false = direct call to
% interp2). If not included, default is OPT.USE_IMPROFILE = false

% Modified Eric M. Moult 07 March 2022

if isempty(varargin)
    OPT.USE_IMPROFILE = false;
else
    OPT = varargin{1};
end

% Take a smaller sub-image around the line.
if (X0(1) <= X0(2))
    X2 = round([X0(1)-1.5*dW, X0(2) + 1.5*dW]);
else
    X2 = round([X0(1)+1.5*dW, X0(2) - 1.5*dW]);
end;

if (Y0(1) <= Y0(2))
    Y2 = round([Y0(1)-1.5*dW, Y0(2) + 1.5*dW]);
else
    Y2 = round([Y0(1)+1.5*dW, Y0(2) - 1.5*dW]);
end;

[ySize, xSize, nChan] = size(img);
    
bdryX = round(max(X2, 1)); bdryX = round(min(bdryX, xSize));  % coordinates of the edges of the ROI
bdryY = round(max(Y2, 1)); bdryY = round(min(bdryY, ySize));
subImg2 = img(min(bdryY):max(bdryY), min(bdryX):max(bdryX), :);

qM = atan2d(diff(Y0), diff(X0));  % angle of the line
perpV = [-sind(qM), cosd(qM)];  % perpendicular unit-vector
L = round((diff(X0)^2 + diff(Y0)^2)^.5);  % length of the line

fiberImg = zeros(2*dW+1, round(L), nChan);

if OPT.USE_IMPROFILE
    for m = -dW:dW;  % Take a series of parallel line-profiles to either side of the fiber.
        X3 = X0 + m*perpV(1) - min(bdryX);
        Y3 = Y0 + m*perpV(2) - min(bdryY);
        fiberImg(m+dW+1,:,:) = squeeze(improfile(subImg2, X3, Y3, L, 'bilinear'));
    end;
else
    interp_points_X = zeros(1, size(fiberImg,1)*size(fiberImg,2));
    interp_points_Y = zeros(1, size(fiberImg,1)*size(fiberImg,2));
    offsets = -dW:dW;
    num_offsets = length(offsets);
    for offset_ind = 1 : num_offsets % Take a series of parallel line-profiles to either side of the fiber.
        X3 = X0 + offsets(offset_ind)*perpV(1) - min(bdryX);
        Y3 = Y0 + offsets(offset_ind)*perpV(2) - min(bdryY);
        interp_points_X((offset_ind - 1)*round(L) + 1:  offset_ind*round(L)) = linspace(X3(1), X3(2), L)';
        interp_points_Y((offset_ind - 1)*round(L) + 1:  offset_ind*round(L)) = linspace(Y3(1), Y3(2), L)';
    end;
    [XX, YY] = meshgrid(1:size(subImg2,2), 1:size(subImg2,1));
    for chan_ind = 1 : nChan
        Vq = interp2(XX,YY,subImg2(:,:,chan_ind), interp_points_X, interp_points_Y, 'bilinear', nan);
        Vq = reshape(Vq,round(L), 2*dW+1)';
        fiberImg(:, :, chan_ind) = Vq;
    end
end
% imagesc(fiberImg); figure; imagesc(Vq)
