function intens = apply_clickydF(roi, movie_in, disp_img,range)
% function intens = apply_clicky(roi, movie_in, disp_img,range)
% Applies an ROI to a movie_in and returns the intensities at those
% locations.  To not display the image, but just return the intensities,
% set the disp_img variable to 'No'.  Default disp_img is yes.
%
% Created: 10/16/11 by AEC and JMK
% 11/4/11: JMK added option to display images

if nargin == 3;
    disp_img = 'yes';
else
    disp_img = lower(disp_img);
end

switch disp_img
    case 'yes'
        nframes = size(movie_in, 3);
        refimg = mean(movie_in, 3);
        
        figure
        subplot(1,3,1)
        imshow(refimg, [], 'InitialMagnification', 'fit')
        hold on;
        
        [ysize, xsize] = size(refimg);
        nroi = length(roi);
        intens = zeros(nframes, nroi);
        
        
        
        colorindex = 0;
        order = get(gca,'ColorOrder');
        [x, y] = meshgrid(1:xsize, 1:ysize);
        for j = 1:nroi;
            xv = roi{j}(:,1);
            yv = roi{j}(:,2);
            inpoly = inpolygon(x,y,xv,yv);
            
            subplot(1,3,1)
            %draw the bounding polygons and label them
            currcolor = order(1+mod(colorindex,size(order,1)),:);
            plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
            text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
            
            itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
            itrace = itrace/mean(itrace(range));
            subplot(1,3,2) % plot the trace
            hold on;
            plot(itrace,'Color',currcolor);
            colorindex = colorindex+1;
            
          subplot(1,3,3) % plot the trace
            hold on;
            plot(mat2gray(itrace),'Color',currcolor);
            colorindex = colorindex+1;
            
            intens(:,j) = itrace;
        end
        hold off
        
    case 'no'
        [ysize, xsize, nframes] = size(movie_in);
        nroi = length(roi);
        intens = zeros(nframes, nroi);
        [x, y] = meshgrid(1:xsize, 1:ysize);

        for j = 1:nroi;
            xv = roi{j}(:,1);
            yv = roi{j}(:,2);
            inpoly = inpolygon(x,y,xv,yv);
            
            itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
            
            intens(:,j) = itrace;
        end
end





