function [roi_points, pixlist1, intens, Vout, corrimg, weightimg, offsetimg] = extractClicky(movie_in, refimg, cmin, cmax);
% function [roi_points, pixlist, intens, Vout, corrimg, weightimg, offsetimg] = extractClicky(movie_in, refimg, cmin, cmax);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% cmin, cmax: contrast range for displaying refimg
% 
% roi_points: Cell array of boundaries of ROIs
% pixlist: list of pixel indices in each ROI
% intens: intensity trace for each ROI
% Vout: extracted Vout trace for each ROI
% corrimg: correlation image (only in pixels specified in pixlist) for each
% ROI
% weightimg: weight image only in pixels specified in pixlist) for each
% ROI
% offsetimg: offset image only in pixels specified in pixlist) for each
% ROI
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Uses extractV to refine the ROI.
% Output Vout is a weighted sum of pixel values, with background
% subtracted.
% Background is defined as the value of the pixel in the ROI with the
% smallest time-average intensity.
% Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% AEC 9 Nov. 2013.

if nargin == 1;
    refimg = mean(movie_in, 3);
end;

if nargin == 2;
    cmin = double(min(refimg(:)));
    cmax = double(max(refimg(:)));
end;

nframes = size(movie_in, 3);
colorimg = repmat((double(refimg)-cmin)/(cmax - cmin), [1 1 3]);
figure
subplot(3,1,1)
imshow(colorimg, 'InitialMagnification', 'fit')
hold on;


[ysize, xsize] = size(refimg(:,:));
npts = 1;
nroi = 1;
intens = [];
Vout = [];
  [x, y] = meshgrid(1:xsize, 1:ysize);
 while(npts > 0)

  subplot(3,1,1)
  [xv, yv] = (getline(gca, 'closed'));
  if size(xv,1) < 3  % exit loop if only a line is drawn
  break
  end
  inpoly = (inpolygon(x,y,xv,yv));
  pixlist1 = find(inpoly(:));
  npix = length(pixlist1);
  pixlistAll = repmat(pixlist1, [1, nframes]) + xsize*ysize*ones(npix, 1)*(0:(nframes-1));
  submov = double(reshape(movie_in(pixlistAll(:)), [npix, nframes]));
  itrace = mean(submov, 1);
  [Vtraj, corrimg1, weightimg1, offsetimg1] = extractV2d(submov, itrace);
  Vtraj = submov'*weightimg1/npix;   % calculate intensity as a weighted mean
  Vtraj = Vtraj - min(mean(submov, 2));  % subtract off the background.
  
  weightmask = zeros(ysize, xsize);
  weightmask(pixlist1) = weightimg1;
  weightmask = mat2gray(double(weightmask));
 
  col = rand(1,3);
  col = col/(sqrt(sum(col.^2)));
  for j = 1:3;
      colorimg(:,:,j) = colorimg(:,:,j) + 0.4*col(j)*(weightmask-.05);
  end;
  subplot(3,1,1);
  imshow(colorimg)

 %draw the bounding polygons and label them
%   plot(xv, yv, 'Linewidth', 1,'Color',col);
  text(mean(xv),mean(yv),num2str(nroi),'Color',col,'FontSize',12);
  
  subplot(3,1,2:3) % plot the trace
  hold on;
  plot(itrace,'Color',col);
  
  intens = [intens; itrace];
  Vout = [Vout Vtraj];
  roi_points{nroi} = [xv, yv];
  pixlist{nroi} = pixlist1;
  corrimg{nroi} = corrimg1;
  weightimg{nroi} = weightimg1;
  offsetimg{nroi} = offsetimg1;
  nroi = nroi + 1;
 end
 intens = intens';
 Vout = Vout;