function [roi_points, intens] = clicky(movie_in, refimg);
% function [roi_points, intens] = clicky(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% AEC and JMK, 12 Oct. 2011.

if nargin < 2;
    refimg = mean(movie_in, 3);
end;

nframes = size(movie_in, 3);

figure
subplot(1,3,1)
imshow(refimg, [], 'InitialMagnification', 'fit') %scales the reference image to fit the clicky image window
hold on;

[ysize, xsize] = size(refimg(:,:,1)); %because it will return the number of rows and then columns
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
  [x, y] = meshgrid(1:xsize, 1:ysize);
 while(npts > 0) %where does npts go to <= 0?

  subplot(1,3,1)
  [xv, yv] = (getline(gca, 'closed')); %getline(fig) lets you select a polyline in the current axes of figure fig using the mouse. 
  %Coordinates of the polyline are returned in X and Y. Use normal button clicks to add points to the polyline. A shift-, right-, or double-click adds a final point and ends the polyline selection. 
  %Pressing Return or Enter ends the polyline selection without adding a final point. Pressing Backspace or Delete removes the previously selected point from the polyline.
  if size(xv,1) < 3  % exit loop if only a line is drawn
  break
  end
  inpoly = inpolygon(x,y,xv,yv);
  
 %draw the bounding polygons and label them
  currcolor = order(1+mod(colorindex,size(order,1)),:); %wrap the colors so that if there are more plots than "order" holds, 
  %Then the current color will wrap around to a color already used
  plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
  text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
  
  itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:)); %repmat(A,[m n p...]) produces a multidimensional array B composed of copies of A. The size of B is [size(A,1)*m, size(A,2)*n, size(A,3)*p, ...].

  subplot(1,3,2:3) % plot the trace
  hold on;
  plot(itrace,'Color',currcolor);
  colorindex = colorindex+1;
  
  intens = [intens; itrace'];
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
 end
 intens = intens';