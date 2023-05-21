function [roi_points, intens] = clickyadjust(movie_in, refimg,total_int_flag);
% function [roi_points, intens] = clicky(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
%
% total_int_flag = 1 will return the total intensity of the highlighted
% region rather than the mean intensity. Default is 0. xoxo sami & simon
%
% AEC and JMK, 12 Oct. 2011.
% BT, 30 Nov. 2021

if nargin < 2 || isempty(refimg)
    refimg = mean(movie_in, 3);
end

if nargin < 3;
    total_int_flag = 0;
end

nframes = size(movie_in, 3);

figure
subplot(1,3,1)
imshow(refimg, [], 'InitialMagnification', 'fit')
hold on;

[ysize, xsize] = size(refimg(:,:,1)); 
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
  [x, y] = meshgrid(1:xsize, 1:ysize);
 while(npts > 0)

  subplot(1,3,1)
  [xv, yv] = (getline(gca, 'closed'));
  if size(xv,1) < 3  % exit loop if only a line is drawn
  break
  end
  inpoly = inpolygon(x,y,xv,yv);
  
 %draw the bounding polygons and label them
  currcolor = order(1+mod(colorindex,size(order,1)),:);
  plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
  text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
  
  if ~total_int_flag
      itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
  else
      itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))));
  end
  
  subplot(1,3,2:3) % plot the trace
  hold on;
  plot(itrace,'Color',currcolor);
  colorindex = colorindex+1;
  xlabel('Frame In Movie')
  ylabel('Intensity')
  title('Intensity vs Frames')
  
  intens = [intens; itrace'];
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
 end
 intens = intens';