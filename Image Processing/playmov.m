function playmov(mov, varargin)
%   Opens a new figure and plays a 3D matrix as a movie; display parameters
%   are set by the variable arguments. The whole matrix is scaled to the
%   maximum and minimum value.
%  
%     playmov(mov, dT) plays the movie, showing a new frame every dT
%     seconds. If dT is not defined, the default is 0.01 seconds.
%  
%     playmov(mov, dT, repeat) plays the movie on loop if repeat = 1. 
%     If repeat is any other value or undefined the movies plays once.
%  
%     playmov(mov, dT, repeat, trange) plays only the frames of the movie
%     specified by trange. trange should be an array. For example, if you 
%     want to play the 10th throug the 20th frame trange = [10:1:20]
%
%     playmov(mov, dT, repeat, trange, colorrange) allows you to set the
%     maximum and minimum values displayed for the movie. colorrange should
%     be a two member array in the format [min max], where both values are 
%     between 0 and 1. The default is the full range, [0 1].
%   
%     2015-08-04, SLF



if nargin < 5
    colorrange = [0 1];
else
    colorrange = varargin{4};
end

if nargin < 4 
    trange = [1:size(mov,3)];
else 
    trange = varargin{3};    
end

if nargin < 3
    repeat = 0;
else
    repeat = varargin{2};
    if ~ismember(repeat, [0 1])
        repeat = 0;
    end    

end

if nargin < 2
    dT = 0.01;
else
    dT = varargin{1};
end

mov_mg = mat2gray(mov);


first = 1;
figure;
while(repeat|first)
    for i = trange;
        imshow2(mov_mg(:,:,i),colorrange)
        title(['Frame ' num2str(i)]);
        pause(dT)
       
    end
    first = 0;
end
