function intens = apply_clicky_faster(roi, movie_in, disp_img, every_n)
% function intens = apply_clicky(roi, movie_in, disp_img)
% Applies an ROI to a movie_in and returns the intensities at those
% locations.  To display the average image and to plot the intensity trace, set disp_image to 'yes'.  
% To not display the image, but just return the intensities,
% set the disp_img variable to 'No'.  Default disp_img is yes.
%
% Created: 10/16/11 by AEC and JMK
% 11/4/11: JMK added option to display images

% Modified:
%   Vicente Parot 2016
%   Cohen Lab - Harvard University

if isa(roi,'vm') % swap order of 2 first args
    temp = roi;
    roi = movie_in;
    movie_in = temp;
end
if isa(movie_in,'vm') % in either case, convert to matrix
    movie_in = movie_in.toimg.data;
end

if ~exist('every_n','var'), every_n = size(movie_in, 3); end
        nframes = size(movie_in, 3);
        refimg = mean(movie_in, 3);
        [ysize, xsize] = size(refimg);

movie_in = tovec(movie_in);

if nargin == 2;
    disp_img = 'yes';
else
    disp_img = lower(disp_img);
end

[~,nroi] = size(roi);

switch disp_img
    case 'yes'

        clf
        subplot(1,3,1)
        imshow(refimg, [], 'InitialMagnification', 'fit')
        hold on;
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
            
%           itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
%           itrace = tovec(movie_in)'*inpoly(:)/nnz(inpoly); % 2016 VP
%           itrace = mean(movie_in(inpoly,:))'; % 2017 VP has bug for 1-pixel ROIs
            itrace = mean(movie_in(inpoly,:),1)'; % 2019 VP
            
            subplot(1,3,2:3) % plot the trace
            hold on;
            plot(reshape(1:nframes,every_n,[]),reshape(itrace,every_n,[]),'Color',currcolor);
            colorindex = colorindex+1;
            
            intens(:,j) = itrace;
        end
        hold off
        
    case 'no'
        intens = zeros(nframes, nroi);
        [x, y] = meshgrid(1:xsize, 1:ysize);

        for j = 1:nroi;
            xv = roi{j}(:,1);
            yv = roi{j}(:,2);
            inpoly = inpolygon(x,y,xv,yv);
            
%           itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
%           itrace = tovec(movie_in)'*inpoly(:)/nnz(inpoly); % 2016 VP
%           itrace = mean(movie_in(inpoly,:))'; % 2017 VP has bug for 1-pixel ROIs
            itrace = mean(movie_in(inpoly,:),1)'; % 2019 VP
            
            intens(:,j) = itrace;
        end
end





