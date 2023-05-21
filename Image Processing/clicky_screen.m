function [roi_points, intens] = clicky_screen(movie_in, refimg)
% function [roi_points, intens] = clicky_screen(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% AEC and JMK, 12 Oct. 2011, KAW Oct. 2013.

if nargin < 2;
    refimg = mean(movie_in, 3);
end;

maxcolors = [0,169,255;
             0,255,0]/255;

nframes = size(movie_in, 3);

figure(1); clf;
subplot(2,1,1)
imshow(refimg, [], 'InitialMagnification', 'fit')
title('click twice: background first, patched cell second');
hold on;

[ysize, xsize] = size(refimg(:,:,1));
colorindex = 0;
order = get(gca,'ColorOrder');
intens = zeros(nframes,2);
roi_points = cell(1,2);
[x, y] = meshgrid(1:xsize, 1:ysize);
for k = 1:2
    subplot(2,1,1)
    [xv, yv] = (getline(gca, 'closed'));
    inpoly = inpolygon(x,y,xv,yv);
    
    % draw the bounding polygons and label them
    currcolor = order(1+mod(colorindex,size(order,1)),:);
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    
    itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
    maxtrace = squeeze(max(max(movie_in.*repmat(inpoly, [1, 1, nframes]))));
    
    subplot(2,1,2) % plot the trace
    hold on;
    plot(itrace,'Color',currcolor);
    plot(maxtrace,'Color',maxcolors(colorindex+1,:));
    legend('avg bkgd','max bkgd','avg cell','max cell')
    colorindex = colorindex+1;
    
    intens(:,k) = itrace';
    roi_points{k} = [xv, yv];
end