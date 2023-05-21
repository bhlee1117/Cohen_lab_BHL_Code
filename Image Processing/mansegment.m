function trackedcells = mansegment(inputcells, h);
% function trackedcells = mansegment(inputcells, h);
% manually add cells to an image to extract their intensities.  Requires a
% cell array, inputcells, to which the new cells are added.  Each element
% of inputcells is a vector of arbitrary length, specifying [xpts ypts]
% defining the boundary of a cell.
%
% Inputcells can be empty, {}.
% h is the handle to the figure showing the cells to be tracked.  The
% figure must be created outside this function

% AEC 17 Feb. 2010

figure(h)
ncells = length(inputcells);
if ncells > 0;
    hold on;
    for j = 1:ncells;
        xpts = inputcells{j}(1,:);
        ypts = inputcells{j}(2,:);
        plot(xpts, ypts, 'r-');
    end;
end;


% Segment the images by hand
xcoords = ones(ysizeT,1)*(1:xsizeT);
ycoords = (1:ysizeT)'*ones(1,xsizeT);

trackedcells = zeros(ysizeT, xsizeT, nexp);

for m = 1;
    colorimg = zeros(ysizeT, xsizeT, 3);
    colorimg(:,:,1) = 2*mat2gray(avgimgsRtrans(:,:,m));
    colorimg(:,:,2) = 2*mat2gray(avgimgsGtrans(:,:,m));
    colorimg(:,:,3) = trackedcells(:,:,m);
%     for j = 1:3;
%         colorimg2(:,:,j) = mat2gray(avgimgsWL(yTrunc,xTrunc,m));
%     end;
    figure(m)
%     imshow([colorimg colorimg2])
    imshow(colorimg, 'InitialMagnification', 'fit')
    k = 40;
    while 1
        figure(m)
        [xpts, ypts] = getline(gca, 'closed');
        in = inpolygon(xcoords, ycoords, xpts, ypts);
        trackedcells(:,:,m) = trackedcells(:,:,m) + in;
        colorimg(:,:,3) = trackedcells(:,:,m);
        imshow(colorimg,  'InitialMagnification', 'fit')
        cell{m,k} = [xpts, ypts];
        k = k + 1;
    end;
end;

ncells = k - 1;
