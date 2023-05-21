function [roi_points, intens] = clicky_rects(movie_in, refimg, separatefigures);
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

if ~exist('refimg','var');
    refimg = mean(movie_in, 3);
end;
% tmi = tovec(tmi);
if ~exist('separatefigures','var');
    separatefigures = true;
end


% nframes = size(movie_in, 3);

if separatefigures
    figure(200); clf
    figure(199); clf
else
    clf
    subplot 131
end
imshow(refimg, [], 'InitialMagnification', 'fit')
hold on;

[ysize, xsize] = size(refimg(:,:,1)); 
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
roi_points = {};
  [x, y] = meshgrid(1:xsize, 1:ysize);
 while(npts > 0)

    if separatefigures
        figure(199)
    else
        subplot(1,3,1)
    end
%   [xv, yv] = (getline(gca, 'closed'));
%   if size(xv,1) < 3  % exit loop if only a line is drawn
%   break
%   end
    if false

        % ellipse instead of rect
        e = imellipse;
        e = e.getVertices;
        xv = e(:,1);
        yv = e(:,2);
    else
      rectpos = getrect;
      xv = rectpos(ones(5,1),1) + [0 rectpos([3 3]) 0 0]';
      yv = rectpos(ones(5,1),2) + [0 0 rectpos([4 4]) 0]';
    end
  if size(unique([xv, yv],'rows'),1) == 1, break, end
%   inpoly = inpolygon(x,y,xv,yv);
  
 %draw the bounding polygons and label them
  currcolor = order(1+mod(colorindex,size(order,1)),:);
  plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
  text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
  
%   itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
%   itrace = tovec(movie_in)'*inpoly(:)/nnz(inpoly);
%     itrace = mean(tmi(inpoly,:))';
%     [min(roipos) max(roipos)-min(roipos)]
    itrace = extract_single_component(movie_in(...
            max(1,floor(rectpos(2))):min(ceil(rectpos(2)+rectpos(4)),size(movie_in,1)),...
            max(1,floor(rectpos(1))):min(ceil(rectpos(1)+rectpos(3)),size(movie_in,2)),:));
    if separatefigures
        figure(200)
    else
        subplot(1,3,2:3) % plot the trace
    end
  hold on;
  plot(itrace,'Color',currcolor);
  colorindex = colorindex+1;
  
  intens = [intens; itrace'];
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
 end
 intens = intens';