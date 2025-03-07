function [BW,maskedImage] = inv_segmentImage(X,vac_ratio)
%segmentImage Segment image using auto-generated code from imageSegmenter app
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the imageSegmenter app. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 18-Jun-2018
%----------------------------------------------------
XX=2^8-1-X;
[his bin]=imhist(XX);
g=1;
while bin(2,1)>10
[his bin]=imhist(XX,g*100);
g=g+1;
end
gg=1;
cum_his=cumsum(his);
while cum_his(gg,1)<vac_ratio*size(XX,1)*size(XX,2)
    gg=gg+1;
    threshold=bin(gg,1);
end


% Create empty mask.
BW = false(size(X,1),size(X,2));
% 영상 이진화 - 수동 임계값
BW = X > 2^8-1-threshold;

% 능동 윤곽선
iterations = 10;
BW = activecontour(X, BW, iterations, 'Chan-Vese');

% disk을(를) 사용하여 마스크 침식
radius = 3;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imerode(BW, se);

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;
end

