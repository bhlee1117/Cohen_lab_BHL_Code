function intens = apply_clicky(roi, movie_in);
% function intens = apply_clicky(roi, movie_in);

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

  subplot(1,3,2:3) % plot the trace
  hold on;
  plot(itrace,'Color',currcolor);
  colorindex = colorindex+1;
  
  intens(:,j) = itrace;
end
hold off



