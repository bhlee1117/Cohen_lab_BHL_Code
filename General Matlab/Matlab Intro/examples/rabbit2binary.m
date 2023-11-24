% convert color image to binary image

R = imread('X:\\Lab\Computer Code\General Matlab\Matlab Intro\examples\Bunny-Wallpapers-bunny-rabbits-128637_1024_768.jpg');
% Red = double(R(:,:,3));
S = std(double(R),0,3);     % calculate the standard deviation.  Only colored pixels will have a large std.


redLevel = 0.5*(max(max(R(:,:,1))));
blueLevel = 0.8*(max(max(R(:,:,3))));
blueMat = double(R(:,:,3));
redMat = double(R(:,:,1));

figure(5);
subplot(1,2,1);
imshow(blueMat,[]);
title('blue mat');
subplot(1,2,2);
imshow(redMat,[]);
title('red mat');
%a means more than an amount of Red
a = redMat>redLevel;
%b means more than an amount of blue
b = blueMat<blueLevel;
% minimumRedLevel = zeros(size(temp,1),size(temp,2));
onlyBunny = redMat.*a.*b;
% thresh = graythresh(S);    % automatically calculate the threshold for converting to binary image
% bw = im2bw(S,thresh);      % generate the binary image

%% 

figure (6);
subplot (1,4,1);
imshow(onlyBunny);
title ('before erosion');


thresh = 0.4*max(onlyBunny(:));
bw = onlyBunny > thresh;
bw = imerode(bw,strel('disk',9));

figure (6);
subplot (1,4,2);
imshow(bw);
title ('before erosion');

thresh = 0.4*max(onlyBunny(:));
bw = onlyBunny > thresh;
bw = imopen(bw,strel('disk',9));

figure (6);
subplot (1,4,3);
imshow(bw);
title ('imopen');

thresh = 0.4*max(onlyBunny(:));
bw = onlyBunny > thresh;
bw = imdilate(bw,strel('disk',9));

figure (6);
subplot (1,4,4);
imshow(bw);
title ('imdilate');

%% 

figure(3)
subplot(2,3,1)
imshow(R)
subplot(2,3,2)
imshow(rgb2gray(R))
title('grayscale')
subplot(2,3,3)
imagesc(S); axis image; colorbar;
title('standard deviation')
subplot(2,3,4)
imagesc(bw); axis image; colorbar;
title('binary image')
subplot(2,3,5)
imshow(onlyBunny,[]); axis image; colorbar;
saveas(gcf,'rabbit_processing.fig');