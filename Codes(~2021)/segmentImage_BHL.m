function [BW,maskedImage] = segmentImage_BHL(X,vac_ratio)
%segmentImage Segment image using auto-generated code from imageSegmenter app
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the imageSegmenter app. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 18-Jun-2018
%----------------------------------------------------


% �ڵ� ����ȭ

[his bin]=imhist(X);
g=1;
while bin(2,1)>10
[his bin]=imhist(X,g*100);
g=g+1;
end
gg=1;
cum_his=cumsum(his);
while cum_his(gg,1)<vac_ratio*size(X,1)*size(X,2)
    gg=gg+1;
    threshold=bin(gg,1);
end

% ���� ����ȭ - ���� �Ӱ谪
BW = X > threshold;

% �ɵ� ������
iterations = 100;
BW = activecontour(X, BW, iterations, 'Chan-Vese');

% disk��(��) ����Ͽ� ����ũ ��â
radius = 1;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imdilate(BW, se);

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;
end
