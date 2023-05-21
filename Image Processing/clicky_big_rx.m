%%%%%%%%%%%%%%%%%%
% Copyright Q-State Biosciences 
%%%%%%%%%%%%%%%%%%
function [roi_points, intens] = clicky_big(movie_in, caxis, refimg);
% function [roi_points, intens] = clicky(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% AEC and JMK, 12 Oct. 2011.

% Modified by Rosalind Xu to allow caxis specs
figure(199); clf
figure(200); clf
if isa(movie_in,'vm')
    movie_in = movie_in.toimg.data;
end

if nargin < 2
    caxis = [];
elseif nargin < 3
    refimg = mean(movie_in, 3);
end;

% nframes = size(movie_in, 3);
movie_in = tovec(movie_in);
figure(199)
imshow(refimg, caxis, 'InitialMagnification', 'fit')
hold on;

[ysize, xsize] = size(refimg(:,:,1)); 
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
  [x, y] = meshgrid(1:xsize, 1:ysize);
 while(npts > 0)

  figure(199)
  [xv, yv] = (getline(gca, 'closed'));
  if size(xv,1) < 3  % exit loop if only a line is drawn
  break
  end
  inpoly = inpolygon(x,y,xv,yv);
  
 %draw the bounding polygons and label them
  currcolor = order(1+mod(colorindex,size(order,1)),:);
  plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
  text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
  
%   itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
% itrace = (inpoly(:)'*tovec(movie_in))';
    itrace = mean(movie_in(inpoly,:))';

  figure(200) % plot the trace
  hold on;
  plot(itrace,'Color',currcolor);
  colorindex = colorindex+1;
  
  intens = [intens; itrace'];
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
 end
 intens = intens';