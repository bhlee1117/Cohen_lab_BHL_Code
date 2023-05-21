function [roi_points, intens1,intens2] = clicky(movie_in1,movie_in2, framerate,refimg);
% function [roi_points, intens] = clicky(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% AEC and JMK, 12 Oct. 2011.

if nargin < 4;
    refimg = mean(movie_in1, 3);
end;


nframes = size(movie_in1, 3);

RefIm=mean(movie_in1,3);
figure(1);
imshow(RefIm,[]);


figure('position',[1, 1, 1080, 1080])
figure(2);
subplot(3,3,1)
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

  subplot(3,3,1)
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

 [intensN1, pbleach] = rem_pbleach(itrace1 , 200);
 [intensN2, pbleach] = rem_pbleach(itrace2 , 200);
 
 

intensN1=medfilt1(intensN1,4);
intensN2=medfilt1(intensN2,4);
intensN1=intensN1(10:end);
intensN2=intensN2(10:end);
  subplot(3,3,i) % plot the trace
  
  t=0:1/framerate:(size(itrace1)-1)/framerate;
  plot(t,itrace1,'b');
  hold on;
plot(t,itrace2,'r');
  title(strcat('Cell', num2str(i-1)));  
  xlabel('Time (s)');
  ylabel('F'); i=i+1;
  colorindex = colorindex+1;
  
  intens1 = [intens1; itrace1];
  intens2 = [intens2; itrace2];
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
 end