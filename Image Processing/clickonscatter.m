function out = clickonscatter(ch1, ch2);
% out = clickonscatter(ch1, ch2); 
% ch1 red, ch2 green!
%  Plots each pixel in scatterplot (intensity in red channel vs. intensity in green channel)
%  Lets you click on different chunks of pixels in the scatterplot and highlights
%  those pixels on the original image
%  -VV, 11/12/12

[ysize xsize nframes] = size(ch1);
if nframes>1
    ch1 = squeeze(mean(ch1,3));
    ch2 = squeeze(mean(ch2,3));
end

%ch1=mat2gray(ch1);
%ch2=mat2gray(ch2);

figure(1)
clf;
while 1;
    subplot(1,2,1)
    plot(ch1(:),ch2(:), 'b.')
    drawnow
    [xv yv] = getline(gca, 'closed');
    hold off
    if length(xv) < 3;
        break
    end;
    plot(ch1(:),ch2(:), 'b.')
    xlabel('Ch1')
    ylabel('Ch2')
    hold all
    plot(xv, yv, 'r-')
    inpoly = inpolygon(ch1(:),ch2(:),xv,yv);
    selectedpix = reshape(inpoly, ysize, xsize);
    subplot(1,2,2)
    colorimg = zeros(ysize, xsize, 3);
    colorimg(:,:,1) = .6*mat2gray(ch1);
    colorimg(:,:,2) = .6*mat2gray(ch2);
    colorimg(:,:,3) = colorimg(:,:,3) + 1.5*selectedpix;
    imshow(colorimg, [], 'InitialMagnification', 'fit')
    
end;

end