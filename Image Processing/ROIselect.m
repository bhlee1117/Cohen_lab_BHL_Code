function [roi_points] = ROIselect(movie_in, refimg);
% function [roi_points] = ROIselect(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
avgimg = mean(movie_in, 3);

if nargin < 2;
    subplots = 1;
    refimg = mean(movie_in, 3);
else
    subplots = 2;
end;
    
figure
subplot(subplots,1,1)
imshow(refimg, [], 'InitialMagnification', 'fit')
subplot(subplots,1,subplots)
imshow(avgimg, [], 'InitialMagnification', 'fit')
hold on;

[ysize, xsize] = size(avgimg(:,:,1));

npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
[x, y] = meshgrid(1:xsize, 1:ysize);
 while(npts > 0)
 [xv, yv] = (getline(gca, 'closed'));
  if size(xv,1) < 3  % exit loop if only a line is drawn
  break
  end
  inpoly = inpolygon(x,y,xv,yv);
 %draw the bounding polygons and label them
  currcolor = order(1+mod(colorindex,size(order,1)),:);
  plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
  text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
 
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
  colorindex = colorindex+1;
 end
