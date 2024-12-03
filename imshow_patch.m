function imshow_patch(image,cax)


figure; clf;
if nargin<2
    cax=[];
end
if iscell(image)
    for z=1:length(image)
nexttile([1 1])
imshow2(image{z},cax)
title(num2str(z))
    end
else
for z=1:size(image,3)
nexttile([1 1])
imshow2(image(:,:,z),cax)
title(num2str(z))
end
end
end