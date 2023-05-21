function makewave_DMD(field, varargin)

% verify user input 
if (mod(size(varargin, 2), 2) ~= 0 || size(varargin, 2) < 2)
    error('makewave_DMD:numEl','Usage: makewave_DMD(fieldofview, object1, timevect1, object2, timevect2 ... etc.)'); 
    return;  
end

% Read in roi objects and make sure all time vectors have length equal to # of frames

for i = 1:(size(varargin,2) / 2)
    roi_vector(i) = varargin{2*i-1};
    if (numel(varargin{2*i}) ~= field.frames)
        error('makewave_DMD:numEl','Usage: time vectors must have same number of elements as the number of frames in the field of view object');  
    else
        time_vect(i,:) = varargin{2*i};
    end    
end

for i = 1:size(roi_vector, 2)
    fmsOfIntrst = find(time_vect(i,:) > 0);
    for j = 1:size(fmsOfIntrst, 2)
        field.add_obj(roi_vector(i), fmsOfIntrst(j));
    end
end


