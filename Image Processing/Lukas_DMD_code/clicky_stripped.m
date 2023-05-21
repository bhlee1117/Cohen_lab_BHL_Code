% Strips some aspects of clicky---notice that roi_points is a cell array
% when it is returned

function roi_points = clicky_stripped(movie_in, refimg)

if nargin < 2;
    refimg = mean(movie_in, 3);
end;

figure(1)
imshow(refimg, [], 'InitialMagnification', 'fit') %scales the reference image to fit the clicky image window
hold on;

npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
 while(npts > 0) 
  [xv, yv] = (getline(gca, 'closed')); 
  if size(xv,1) < 3  % exit loop if only a line is drawn
  break
  end
  
  %draw the bounding polygons and label them
  currcolor = order(1+mod(colorindex,size(order,1)),:); %wrap the colors 
  
  plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
  text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
 
  colorindex = colorindex+1;
  
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
 end