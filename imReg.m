function [RegImg, tformReg]=imReg(Img2Reg, regImg)

points1 = detectSURFFeatures(mat2gray(Img2Reg));
points2 = detectSURFFeatures(mat2gray(regImg));

[features1, valid_points1] = extractFeatures(mat2gray(Img2Reg), points1);
[features2, valid_points2] = extractFeatures(mat2gray(regImg), points2);

indexPairs = matchFeatures(features1, features2);
matchedPoints1 = valid_points1(indexPairs(:, 1));
matchedPoints2 = valid_points2(indexPairs(:, 2));

try
tformReg = estimateGeometricTransform(matchedPoints1, matchedPoints2, 'affine');
RegImg = imwarp(Img2Reg, tformReg, 'OutputView', imref2d(size(regImg)));
imshowpair(image2, RegImg, 'montage');
catch

  figure(2); clf;
nexttile([1 1]);
imshow2(Img2Reg,[]); hold on;
nexttile([1 1]);
imshow2(regImg,[]); hold on;

points = [];
disp('Click to add points. Press Enter to finish.');
while true
    w = waitforbuttonpress;

    % Check if the key pressed is 'Enter' (character 13)
    if w == 1 && get(gcf, 'CurrentCharacter') == 13
        break;
    end
    h=drawpoint(gca);
    points=[points; h.Position];
end
close(figure(2));

figure(2);
nexttile([1 1]);
imshow2(regImg,[]); hold on;
plot(points(:,1),points(:,2),'ro')
text(points(:,1),points(:,2)+10,num2str([1:size(points,1)]'),'color','w');

ax2=nexttile([1 1]);
imshow2(Img2Reg,[]); hold on;
points2Reg = [];
disp('Click to add points. Press Enter to finish.');
while true
    w = waitforbuttonpress;

    % Check if the key pressed is 'Enter' (character 13)
    if w == 1 && get(gcf, 'CurrentCharacter') == 13
        break;
    end
    h=drawpoint(gca);
    points2Reg=[points2Reg; h.Position];
end
close(figure(2));

tformReg = estimateGeometricTransform(points2Reg, points, 'affine');
RegImg = imwarp(Img2Reg, tformReg, 'OutputView', imref2d(size(regImg)));
figure(2); clf;
nexttile([1 1])
imshowpair(regImg, RegImg, 'montage');
nexttile([1 1])
imshow2(imfuse(regImg,RegImg),[])
end

end