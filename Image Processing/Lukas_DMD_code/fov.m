% Create a field of view (fov) and specify functions that add objects to
% the field of view

classdef fov < handle
    properties
        height
        width
        frames
        data
    end
    methods
        
        function fov = fov(h, l, f)
            fov.height = h; 
            fov.width = l; 
            fov.frames = f; 
            fov.data = zeros(h, l, f); 
        end
        
        function add_obj(fov, obj, frame)
            x = fov.data;
        
            object = obj.data;
            
            %Deal with coordinates of the object being out of the field of
            %view
            if ((obj.top_lefty + obj.height - 1) > fov.height)
                yindices = obj.top_lefty:fov.height; 
                object = object(1:numel(yindices), :);
            elseif (obj.top_lefty < 1)
                yindices = 1:(obj.height - abs(obj.top_lefty) - 1); 
                object = object((abs(obj.top_lefty) + 2):end,:); 
            else
                yindices = obj.top_lefty:(obj.top_lefty + obj.height - 1); 
            end   

            
            if ((obj.top_leftx + obj.width - 1) > fov.width)
                xindices = obj.top_leftx:fov.width; 
                object = object(:, 1:numel(xindices)); 
            elseif (obj.top_leftx < 1)
                xindices = 1:(obj.width - abs(obj.top_leftx) - 1); 
                object = object(:,(abs(obj.top_leftx) + 2):end); 
            else
                xindices = obj.top_leftx:(obj.top_leftx + obj.width - 1); 
            end
                
            % Add the object on top of any other objects that may already
            % be in the field of view   
            x(yindices, xindices, frame) = x(yindices, xindices, frame) + object;
            fov.data = (x > 0); 
        end
    end
end