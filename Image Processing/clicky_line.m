function [roi_points, intens] = clicky(movie_in,linewidth,refimg);
% function [roi_points, intens] = clicky(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% linewidth: fine-tune the width of the line ROI
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% He Tian, modified from clicky.m, 02/26/2021

if nargin < 3;
    refimg = mean(movie_in, 3);
end;

if nargin<2
    linewidth=2;
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

while(npts > 0)
      subplot(1,3,1)
      h=drawpolyline(gca);
      xv=round(h.Position(:,1));
      yv=round(h.Position(:,2)); 
      
      mask=zeros(ysize,xsize);
      if length(xv) < 2  % exit loop if only a point is drawn
         h.Visible=0;
         break
      end
      % connect the dots to create a continuous countour
      for i=1:length(xv)-1
          mask(yv(i),xv(i))=1;
      
          if xv(i+1)~=xv(i)
            slope=(yv(i+1)-yv(i))/(xv(i+1)-xv(i));
            intercept=yv(i)-slope*xv(i);
              for j=min(xv(i),xv(i+1))+1:max(xv(i),xv(i+1))
                  mask(round(slope*j+intercept),j)=1;
              end
          else
              for j=min(xv(i),xv(i+1))+1:max(xv(i),xv(i+1))
              mask(yv(i),j)=1;
              end
          end
      end
      %mask=bwmorph(mask,'bridge');
      se=strel('sphere',linewidth);
      mask2=imdilate(mask,se);
% 
%       figure(2);clf;
%       subplot(2,1,1)
%       imshow(mask)
%       subplot(2,1,2)
%       imshow(mask2)

      % get the coordinates of the masked pixels
      pixels=find(mask2==1);
      [row,col] = ind2sub(size(mask2),pixels);
    
    currcolor = order(1+mod(colorindex,size(order,1)),:);
     %draw the bounding polygons and label them
     subplot(1,3,1)
     h.Visible=0;
     plot(col, row, '.','MarkerSize', 1,'Color',currcolor);
     text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12); 
      currcolor = order(1+mod(colorindex,size(order,1)),:); %wrap the colors so that if there are more plots than "order" holds, 
      %Then the current color will wrap around to a color already used

      itrace = squeeze(sum(sum(movie_in.*mask2)))/length(pixels); %repmat(A,[m n p...]) produces a multidimensional array B composed of copies of A. The size of B is [size(A,1)*m, size(A,2)*n, size(A,3)*p, ...].

     
     subplot(1,3,2:3) % plot the trace
     plot(itrace,'Color',currcolor);
     colorindex = colorindex+1;
     hold on
      intens = [intens; itrace'];
      roi_points{nroi} = [col, row];
      nroi = nroi + 1;
end
  
 intens = intens';