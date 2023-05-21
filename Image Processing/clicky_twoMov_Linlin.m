function [roi_points, intens1,intens2] = clicky_twoMov_Linlin(movie_in1,movie_in2, framerate,ref);

if ref == 1;
    refimg = mean(movie_in1, 3);
else
    refimg = mean(movie_in2,3);
end;

nframes = size(movie_in1, 3);

figure;
subplot(3,1,1)
imshow(refimg, [], 'InitialMagnification', 'fit')
hold on;

[ysize, xsize] = size(refimg(:,:,1)); 
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens1 = [];
intens2 = [];
  [x, y] = meshgrid(1:xsize, 1:ysize);
  i=2;
 while(npts > 0)

  subplot(3,1,1)
  [xv, yv] = (getline(gca, 'closed'));
  if size(xv,1) < 3  % exit loop if only a line is drawn
  break
  end
  inpoly = inpolygon(x,y,xv,yv);
  
 %draw the bounding polygons and label them
  currcolor = order(1+mod(colorindex,size(order,1)),:);
  plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
  text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
  
  itrace1 = squeeze(sum(sum(movie_in1.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
  itrace2 = squeeze(sum(sum(movie_in2.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));

 [intensN1, pbleach] = rem_pbleach(itrace1 , 100);
 [intensN2, pbleach] = rem_pbleach(itrace2 , 100);

  subplot(3,1,2) % plot the trace
  
  t=0:1/framerate:(size(itrace1)-1)/framerate;
  plot(t,itrace1);
  title('Andor green');  
  xlabel('Time (s)');
  ylabel('F'); i=i+1;
  hold on;
  
  subplot (3,1,3)
  plot(t,itrace2);
  title('HM red');  
  xlabel('Time (s)');
  ylabel('F');
  hold on;
  
  colorindex = colorindex+1;
  
  intens1 = [intens1; itrace1'];
  intens2 = [intens2; itrace2'];
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
 end