% Interactive Line-Mask Explorer in MATLAB
% Usage: Call interactive_line2mask(image, masks)
% image: s1 x s2 grayscale image
% masks: s1 x s2 x N logical mask array

function sequence=interactive_line2mask(image, masks)
    figure;
    imshow(image, []);
    title('Draw a line: click two points');

    % User draws line
    [x, y] = ginput(2);
    x = round(x);
    y = round(y);

    % Generate linear indices along the line
    [rr, cc] = drawline(y(1), x(1), y(2), x(2));
    line(cc, rr, 'Color', 'g');

    % Track mask index at each point
    mask_sequence = zeros(length(rr), 1);

    for i = 1:length(rr)
        for j = 1:size(masks, 3)
            if masks(rr(i), cc(i), j)
                mask_sequence(i) = j;
                break;  % take the first mask that matches
            end
        end
    end

    % Remove consecutive duplicates to show entry sequence
    sequence = mask_sequence(mask_sequence > 0);
    sequence = sequence([true; diff(sequence) ~= 0]);

    fprintf('Sequence of mask indices along the line: ');
    disp(sequence');
end

function [rr, cc] = drawline(r1, c1, r2, c2)
    % Bresenham's line algorithm
    [rr, cc] = bresenham(r1, c1, r2, c2);
    valid = rr >= 1 & rr <= size(getimage, 1) & cc >= 1 & cc <= size(getimage, 2);
    rr = rr(valid);
    cc = cc(valid);
end

function [x, y] = bresenham(x1, y1, x2, y2)
    x1 = round(x1); y1 = round(y1); x2 = round(x2); y2 = round(y2);
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    steep = dy > dx;

    if steep
        [x1, y1] = deal(y1, x1);
        [x2, y2] = deal(y2, x2);
        [dx, dy] = deal(dy, dx);
    end

    if x1 > x2
        [x1, x2] = deal(x2, x1);
        [y1, y2] = deal(y2, y1);
    end

    derr = 2 * dy;
    err = 0;
    ystep = 1;
    if y1 > y2
        ystep = -1;
    end
    y = y1;
    x_vec = [];
    y_vec = [];

    for x = x1:x2
        if steep
            x_vec(end+1) = y;
            y_vec(end+1) = x;
        else
            x_vec(end+1) = x;
            y_vec(end+1) = y;
        end
        err = err + derr;
        if err > dx
            y = y + ystep;
            err = err - 2*dx;
        end
    end

    x = y_vec(:);
    y = x_vec(:);
end
