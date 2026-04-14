function [RegImg, tformReg] = interactiveImageRegistration(Img2Reg, refImg)
% INTERACTIVEIMAGEREGISTRATION Align Img2Reg to refImg using click + refinement
%   Outputs registered image and transformation

    % Init
    close all;
    RegImg = [];
    tformReg = [];    
    fig = figure('Name','Interactive Image Registration','KeyPressFcn',@keyPress);

    % Axes
    tiledlayout(3,2);
    ax1 = nexttile([2 1]); imshow(Img2Reg,[prctile(Img2Reg(:),1) prctile(Img2Reg(:),99)]); title('Image to Register'); hold on;
    ax2 = nexttile([2 1]); imshow(refImg,[prctile(refImg(:),1) prctile(refImg(:),99)]); title('Reference Image'); hold on;
    ax3 = nexttile([1 2]); hOverlay = imshow(zeros(size(refImg)),[]); title('Registered Overlay'); hold on;

    points1 = [];
    points2 = [];
    txtHandles1 = [];
    txtHandles2 = [];

    % Instructions
    disp('Click matching points: First in Img2Reg (Top), then in RefImg (Middle). Minimum 3 pairs. Press "q" to finish. Press "u" to undo last pair.');

    % Click loop
    while true
        try
            % First click on Img2Reg
            figure(fig); axes(ax1);
            h1 = drawpoint;
            if isempty(h1); break; end
            pt1 = h1.Position;
            %pt1 = refineLocalXCorr(Img2Reg, pt1);
            points1(end+1,:) = pt1;
            delete(h1);
            plot(ax1, pt1(1), pt1(2), 'ro');
            txtHandles1(end+1) = text(ax1, pt1(1)+5, pt1(2), num2str(size(points1,1)), 'Color', 'w', 'FontSize', 10);

            % Then click on RefImg
            axes(ax2);
            h2 = drawpoint;
            if isempty(h2); break; end
            pt2 = h2.Position;
            %pt2 = refineLocalXCorr(refImg, pt2);
            points2(end+1,:) = pt2;
            delete(h2);
            plot(ax2, pt2(1), pt2(2), 'ro');
            txtHandles2(end+1) = text(ax2, pt2(1)+5, pt2(2), num2str(size(points2,1)), 'Color', 'w', 'FontSize', 10);

            % If enough points, register and show overlay
            updateRegistration();
        catch
            break;
        end
    end

    function ptRefined = refineLocalXCorr(im, pt)
        % Look for subpixel match in 11x11 region
        w = 5; % half window size
        x = round(pt(1));
        y = round(pt(2));
        x = min(max(x, w+1), size(im,2)-w);
        y = min(max(y, w+1), size(im,1)-w);
        template = im(y-w:y+w, x-w:x+w);
        level = graythresh(template);
        template_bin=template>level*2;

        c = normxcorr2(template_bin, im);
        [~, maxIdx] = max(c(:));
        [ypeak, xpeak] = ind2sub(size(c), maxIdx);
        ptRefined = [xpeak-w-1, ypeak-w-1];
    end

    function updateRegistration()
        if size(points1,1) >= 3 && size(points2,1) >= 3
            tformReg = estimateGeometricTransform2D(points1, points2, 'affine');
            RegImg = imwarp(Img2Reg, tformReg, 'OutputView', imref2d(size(refImg)));
            axes(ax3);
            imshowpair(mat2gray(refImg), mat2gray(RegImg), 'Scaling', 'joint', 'Parent', ax3);
            title(ax3, 'Registered Overlay');
        end
    end

    function keyPress(~, event)
        switch event.Key
            case 'q'
                uiresume;
                close(fig);
            case 'u'
                if ~isempty(points1)
                    points1(end,:) = [];
                    delete(txtHandles1(end));
                    txtHandles1(end) = [];
                end
                if ~isempty(points2)
                    points2(end,:) = [];
                    delete(txtHandles2(end));
                    txtHandles2(end) = [];
                end
                cla(ax3); % clear overlay
                updateRegistration();
        end
    end
end
