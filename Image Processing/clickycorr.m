function [roi_points, intens] = clickycorr(movie_in, refimg);
% function [roi_points, intens] = clicky(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete. The function will wait 
% for keyboard input. For an autocorrelation of the selected area, input a. 
% For a cross correlation, input x and click another region.
% Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the correlation traces.
% AEC and JMK, 12 Oct. 2011.

if nargin < 2;
    refimg = mean(movie_in, 3);
end;

nframes = size(movie_in, 3);

figure
subplot(1,3,1)
imshow(refimg, [], 'InitialMagnification', 'fit');
title('a for autocorrelation, x for cross correlation');
hold on;

[ysize, xsize] = size(refimg(:,:,1)); 
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
[x, y] = meshgrid(1:xsize, 1:ysize);

xcorrflag = 0;
while(npts > 0)

    subplot(1,3,1)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
        break
    end
    
    if xcorrflag ==0;
        while(1)
            waitforbuttonpress;
            input = get(gcf,'CurrentCharacter');
            if input =='a' || input=='x'
                break
            end
        end
    end
    
    inpoly = inpolygon(x,y,xv,yv);

    %draw the bounding polygons and label them
    if xcorrflag
        currcolor = order(1+mod(colorindex,size(order,1)),:);
        text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    else
        currcolor = order(1+mod(colorindex,size(order,1)),:);
        text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    end
    
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
    
    if input=='a'
        [Corr,lags] = xcorr(itrace-mean(itrace));
        subplot(1,3,2:3) % plot the trace
        hold on;
        plot(lags, Corr,'Color',currcolor);
        colorindex = colorindex+1;

        intens = [intens; Corr'];
        roi_points{nroi} = [xv, yv];
        nroi = nroi + 1;
    elseif input =='x' && xcorrflag==0;
        itraceprev=itrace;
        xcorrflag = 1;
    elseif xcorrflag
        [Corr,lags] = xcorr(itrace-mean(itrace),itraceprev-mean(itraceprev));
        subplot(1,3,2:3) % plot the trace
        hold on;
        plot(lags,Corr,'Color',currcolor);
        colorindex = colorindex+1;

        intens = [intens; Corr'];
        roi_points{nroi} = [xv, yv];
        nroi = nroi + 1;
        xcorrflag = 0;
    end  
end
intens = intens';
 
 
 
 
 