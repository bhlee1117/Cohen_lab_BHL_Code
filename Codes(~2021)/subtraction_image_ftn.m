%  Autofluorescence subtraction code
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code removes the autofluorescent signal from red channel (channel 1)
% from green channel (channel 2) accuracy exceed 90%.

% INPUTS 

% Two channel image matrix size(im,3)=2
% 

% OUTPUTS
% Subtracted image matrix

% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.

function subtraction=subtraction_image_ftn(im)
        stack=double(im);
        R=stack(:,:,1);
        
         if max(max(stack(:,:,1)))<125
        subtraction=stack(:,:,2);
         else
        subtraction=stack(:,:,2)-R/max(max(R))*max(max(stack(:,:,2)));
         end
        subtraction=uint16(subtraction);
end
