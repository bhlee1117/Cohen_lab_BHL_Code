function [roi_points, intens] = nested_clicky_rects(movie_in, refimg);
% function [roi_points, intens] = clicky(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% AEC and JMK, 12 Oct. 2011.

% 2016 Vicente Parot
% Cohen Lab - Harvard University

if isa(movie_in,'vm')
    movie_in = movie_in.toimg.data;
end

if nargin < 2;
refimg = mean(movie_in, 3);
end;

figure(200), clf reset
figure(199), clf reset
imshow(refimg, [], 'InitialMagnification', 'fit')
hold on;
drawnow

[ysize, xsize] = size(refimg(:,:,1)); 
npts = 1;
order = get(gca,'ColorOrder');
intens = [];
roi_points = {};
while(npts > 0)
%     subplot(1,3,1)
    figure(199)
	rectpos = getrect;
    xv = rectpos(ones(5,1),1) + [0 rectpos([3 3]) 0 0]';
    yv = rectpos(ones(5,1),2) + [0 0 rectpos([4 4]) 0]';
    if size(unique([xv, yv],'rows'),1) == 1, break, end

    submov = movie_in(...
        max(1,ceil(rectpos(2))):min(ysize,floor(rectpos(2)+rectpos(4))),...
        max(1,ceil(rectpos(1))):min(xsize,floor(rectpos(1)+rectpos(3))),...
        :);
    subimg = refimg(...
        max(1,ceil(rectpos(2))):min(ysize,floor(rectpos(2)+rectpos(4))),...
        max(1,ceil(rectpos(1))):min(xsize,floor(rectpos(1)+rectpos(3))),...
        :);
    [subroipts, subintens(subintens)] = clicky_rects(submov,subimg);
    for it = 1:numel(subroipts)
        subroipts{it}(:,1) = subroipts{it}(:,1) + ceil(rectpos(1)) - 1;
        subroipts{it}(:,2) = subroipts{it}(:,2) + ceil(rectpos(2)) - 1;
    end
    intens = [intens subintens];
    roi_points = [roi_points subroipts];

%     clf
%     subplot(1,3,1)
    clf(199)
    imshow(refimg, [], 'InitialMagnification', 'fit')
    hold on;

    colorindex = 0;
    for it = 1:numel(roi_points)
%         subplot 131
        figure(199)
        currcolor = order(1+mod(colorindex,size(order,1)),:);
        plot(roi_points{it}(:,1), roi_points{it}(:,2), 'Linewidth', 1,'Color',currcolor);
        text(mean(roi_points{it}(:,1)),mean(roi_points{it}(:,2)),num2str(colorindex+1),'Color',currcolor,'FontSize',12);

%         subplot(1,3,2:3) % plot the trace
        figure(200)
        plot(intens(:,it),'Color',currcolor);
        hold on
        colorindex = colorindex+1;
    end
    drawnow
end
