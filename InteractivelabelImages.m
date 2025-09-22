function selection = InteractivelabelImages(imgCell)
% InteractivelabelImages Interactive labeling of images in a cell array.
%   selection = InteractivelabelImages(imgCell)
%   - imgCell: cell array of images (e.g. grayscale or RGB, various sizes)
%   - selection: binary vector (1 = selected with Enter, 0 = skipped with 's')

n = numel(imgCell);
selection = NaN(1, n);  % preallocate
i = 1;

fig = figure('Name', 'Image Labeling', ...
    'KeyPressFcn', @onKeyPress, ...
    'CloseRequestFcn', @onClose);

updateImage();
uiwait(fig);  % <- wait for user input to complete

    function updateImage()
        if ishandle(fig)
            clf;
            imagesc(imgCell{i});
            %axis image off
            title(sprintf('Image %d / %d\n[Enter]=Mark as 1 | [s]=Skip | [b]=Back', i, n), 'FontSize', 14);
        end
    end

    function onKeyPress(~, event)
        switch event.Key
            case 'return'  % Enter
                selection(i) = 1;
                i = i + 1;
            case 's'  % skip
                selection(i) = 0;
                i = i + 1;
            case 'b'  % back
                i = max(i - 1, 1);
            otherwise
                return
        end

        if i > n
            uiresume(fig);
            close(fig);
        else
            updateImage();
        end
    end

    function onClose(~, ~)
        % allow clean closure
        uiresume(fig);
        delete(fig);
    end
end
