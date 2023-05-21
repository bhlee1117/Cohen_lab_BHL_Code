function out = clickonscatteradjust(ch1, ch2, originim, data1, data2);
% out = clickonscatter(ch1, ch2); 
% ch1 red, ch2 green!
%  Plots each pixel in scatterplot (intensity in red channel vs. intensity in green channel)
%  Lets you click on different chunks of pixels in the scatterplot and highlights
%  those pixels on the original image
%  -VV, 11/12/12
% Now it draws on the image originim that the user wants, and it shows a best fit linear line to look for outliers
% data1 and data2 are used for the best fit line while ch1 and ch2 are used
% for the highlighting part on the original image
% -BT, 7/15/21

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
    
    x = data1';
    y = data2';
    format long
    b1 = x\y;
    yCalc1 = b1*x;
    scatter(x,y,'b.')
    hold on
    plot(x,yCalc1)
    grid on
    hold off
    
%    plot(ch1(:),ch2(:), 'b.')
    drawnow
    [xv yv] = getline(gca, 'closed');
    hold off
    if length(xv) < 3;
        break
    end;
%    plot(ch1(:),ch2(:), 'b.')

    x = data1';
    y = data2';
    format long
    b1 = x\y;
    yCalc1 = b1*x;
    scatter(x,y,'b.')
    hold on
    plot(x,yCalc1)
    grid on
    hold off

    xlabel('Ch1')
    ylabel('Ch2')
    hold all
    plot(xv, yv, 'r-')
    inpoly = inpolygon(ch1(:),ch2(:),xv,yv);
    selectedpix = reshape(inpoly, ysize, xsize);
    

    
%     subplot(1,2,2)
%     colorimg = zeros(ysize, xsize, 3);
%     colorimg(:,:,1) = 0.2*mat2gray(ch1);
%     colorimg(:,:,2) = 0.2*mat2gray(ch2);
%     colorimg(:,:,3) = selectedpix;
%     imshow(colorimg, [], 'InitialMagnification', 'fit')
    
    subplot(1,2,2)
    colorimg2 = zeros(ysize, xsize, 3);
    colorimg2(:,:,1) = originim;
    colorimg2(:,:,2) = originim+selectedpix;
    colorimg2(:,:,3) = originim;
    imshow(colorimg2, [], 'InitialMagnification', 'fit')
%     
end

end