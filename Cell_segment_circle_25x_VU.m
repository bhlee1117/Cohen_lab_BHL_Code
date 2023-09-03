function [centers, radii]=Cell_segment_circle_25x_VU(avgImg,thres)
%% Control Tower
if nargin<4

end

%%

avgImg=mat2gray(avgImg-medfilt2(avgImg,[150 150],"symmetric"));
im_G=imgaussfilt(avgImg,1.5);

%%
detectCircles = @(x) imfindcircles(x,[6 10], ...
    'Sensitivity',thres, ...
    'EdgeThreshold',0.12, ...
    'Method','TwoStage', ...
    'ObjectPolarity','Bright');
[centers, radii, metric] = detectCircles(im_G);
% imagesc(255-imfilt)
% colormap('gray')
% hold all
% for i=1:size(centers,1)
%     plot(centers(i,1),centers(i,2),'ro')
%     text(centers(i,1),centers(i,2)+10,num2str(i),'color','r',...
% 'fontweight','bold')
% end
%
% axis equal
end