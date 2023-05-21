% Create an roi object. An roi object is a rectangle that bounds some
% polygonal shape specified by roi_points. 

classdef roi < handle
    properties 
        top_leftx
        top_lefty
        height
        width
        data
        coords
    end
    methods
        function roi = roi(varargin) %takes roi_points and can also take a transformation matrix
            if nargin > 2
                error('roi:numEl', 'Usage: roi = roi(roi_points, tform'); 
            elseif nargin == 1
                roi.coords = varargin{1};
                roi.coords = roi.coords{1};
                xv = roi.coords(:,1); 
                yv = roi.coords(:,2); 
            else
                roi_points = varargin{1}; 
                tformMat = varargin{2}; 
                
                temp = roi_points{1}; 
                temp = [temp'; ones(1, size(temp,1))]; 
                temp = tformMat*temp;
                
                xv = temp(1,:)'; 
                yv = temp(2,:)'; 
                
                roi.coords = [xv, yv]; 
            end             
            
            % Define the bounding box surrounding the object
            roi.top_leftx = floor(min(xv)); 
            roi.top_lefty = floor(min(yv));
            roi.width = ceil(max(xv)) - roi.top_leftx;
            roi.height = ceil(max(yv)) - roi.top_lefty; 
            
            % Now create an object that contains the polygon
            [x, y] = meshgrid(1:roi.width, 1:roi.height); 
            inpoly = inpolygon(x, y, xv - roi.top_leftx, yv - roi.top_lefty); 
            
            roi.data = inpoly; 
        end
        
        % updates the x coordinates of the roi objects when the top left
        % value changes
        function set.top_leftx(roi, value)
            if (numel(roi.top_leftx) > 0)
                update_xcoords(roi, value); 
            end
            roi.top_leftx = value; 
        end
        
        % updates the y coordinates of the roi objects when the top left
        % value changes
        function set.top_lefty(roi, value)
            if (numel(roi.top_lefty) > 0)
                update_ycoords(roi, value);
            end
            roi.top_lefty = value; 
        end
        
        % internal function that does arithmitic on the x coord update
        function update_xcoords(roi, value)
            roi.coords(:,1) = roi.coords(:,1) + (value - roi.top_leftx); 
        end
        
        
        % internal function that does arithmitic on the y coord update
        function update_ycoords(roi, value)
            roi.coords(:,2) = roi.coords(:,2) + (value - roi.top_lefty); 
        end
        
        % applies the transformation to the roi object
        function apply_tform(roi, tform)
            
            % apply the transform
            temp = [roi.coords'; ones(1, size(roi.coords,1))]; 
            temp = tform*temp;
            
            xv = temp(1,:)';
            yv = temp(2,:)';
            
            % define the new roi coordinates
            roi.coords = [xv, yv]; 
            
            % Define the bounding box surrounding the object
            roi.top_leftx = floor(min(xv)); 
            roi.top_lefty = floor(min(yv));
            roi.width = ceil(max(xv)) - roi.top_leftx;
            roi.height = ceil(max(yv)) - roi.top_lefty; 
            
            % Now create an object that contains the polygon
            [x, y] = meshgrid(1:roi.width, 1:roi.height); 
            inpoly = inpolygon(x, y, xv - roi.top_leftx, yv - roi.top_lefty); 
            
            roi.data = inpoly; 
        end
    end
end
        