function roilabel = interactive_ROIlabel(imageStack, referenceImage)

    % interactiveLabeling allows interactive labeling of ROIs in an X x Y x N binary image stack.
    % The user labels ROIs starting from label 1, iteratively moving to the next label.
    %
    % Inputs:
    %   - imageStack: An X x Y x N binary image stack.
    %   - referenceImage: An X x Y image to be used as a reference for labeling.
    %
    % Output:
    %   - updatedStack: An updated X x Y x N binary image stack with user-labeled ROIs.

    % Validate inputs
    

    [X, Y, N] = size(imageStack);
    if ~isequal(size(referenceImage), [X, Y])
        error('referenceImage dimensions must match the first two dimensions of imageStack.');
    end

    % Initialize updated stack
    updatedStack = imageStack;

    % Create a figure for interactive labeling
    fig=figure;
    label = 1;
    roilabel=NaN(1,size(imageStack,3));
    colors = lines(100); % Predefined set of colors for labels

    while true
        labeledstack(:,:,label)=false(X, Y);
        clf;
        imshow(referenceImage, []); % Show reference image
        hold on;
        % Overlay contours of all masks in the current stack
        for i = 1:N
            mask = updatedStack(:, :, i);
            [contours, ~] = bwboundaries(mask, 'noholes');
            for k = 1:length(contours)
                plot(contours{k}(:, 2), contours{k}(:, 1), 'r', 'LineWidth', 1);
            end
        end

        if label>1
            newColoredMask = im_merge(labeledstack(:, :, 1:label),colors(1:label,:));
                h = imshow(newColoredMask);
                set(h, 'AlphaData', 0.5);
        end
        hold off;

        % Prompt user to select ROIs
        title(sprintf('Label %d: Click points for ROIs, press Enter to finish, or type ''q'' to quit.', label));
        drawnow;

        % Collect ROI mask for the current label
        currentLabelMask = false(X, Y);        
        while true
            [x, y, button] = ginput(1); % Get a single point input

            if isempty(button) % If Enter is pressed, finish this label
                break;
            elseif button == 'q' % If 'q' is pressed, quit
                close(fig);
                return;
            end            

            % Round coordinates and mark the ROI
            x = round(x);
            y = round(y);
            if (x > 0 && x <= Y && y > 0 && y <= X) & length(find(squeeze(imageStack(y,x,:))>0))==1

                roilabel(find(squeeze(imageStack(y,x,:))>0,1))=label;
                currentLabelMask = double(imageStack(:,:,find(squeeze(imageStack(y,x,:))>0,1))>0);
                labeledstack(:, :, label) = labeledstack(:, :, label) | currentLabelMask;

                % Update display
                imshow(referenceImage, []);
                hold on;

                % Overlay existing contours
                for i = 1:N
                    mask = updatedStack(:, :, i);
                    [contours, ~] = bwboundaries(mask, 'noholes');
                    for k = 1:length(contours)
                        plot(contours{k}(:, 2), contours{k}(:, 1), 'r', 'LineWidth', 1);
                    end
                end

                % Show the new ROI filled with the current label color
                newColoredMask = im_merge(labeledstack(:, :, 1:label),colors(1:label,:));
                h = imshow(newColoredMask);
                set(h, 'AlphaData', 0.5);
                hold off;
            else
                fprintf('Point out of bounds. Try again.\n');
            end
        end

        % Increment label index
        label = label + 1;

        % Check if label exceeds stack size
        if label > N
            fprintf('All slices in the stack have been labeled. Exiting.\n');
            break;
        end
    end

end
