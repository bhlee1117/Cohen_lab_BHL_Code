function roi_points = click_im(im,clicky_type,fig)
if ~exist('fig','var'),fig=figure; end
figure(fig);clf
imshow(im,[],'initialmagnification','fit')
switch clicky_type
    case 'clicky'
        
        % start cliky loop
        hold on
        title('Please clicky ROIs')
        npts = 1;
        nroi = 1;
        while(npts > 0)

            [xv, yv] = (getline(gca, 'closed'));
            if size(xv,1) < 3  % exit loop if only a line is drawn
            break
            end

            %draw the bounding polygons and label them
            plot(xv, yv, 'Linewidth', 1);

            roi_points{nroi} = [xv, yv];
            nroi = nroi + 1;
        end

        
        
    case 'rect'
        
        try %#okay uncaught
            myrect = getrect;   
        end
        if exist('myrect','var') && isnumeric(myrect) && numel(myrect)==4
            roi_points = {myrect(1:2)+myrect(3:4).*[0 0 1 1 0; 0 1 1 0 0]'};
            hold on
            plot(roi_points{1}(:,1),roi_points{1}(:,2))
            hold off
        end
        
    case 'pts'
        [xv, yv] = (getpts(gca, 'closed'));
        roi_points = [xv yv];
    case 'Circle'
    case 'Load Picture'
end
