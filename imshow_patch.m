function imshow_patch(image,cax)


clf;

if iscell(image)
    for z=1:length(image)

        if nargin<2
            cax=[prctile(image{z}(:),1) prctile(image{z}(:),99)];
        end
        nexttile([1 1])
        try
            imshow2(image{z},cax)
        catch
            imshow2(image{z})
        end
        title(num2str(z))
    end
else
    for z=1:size(image,3)

        if nargin<2
            cax=[prctile(tovec(image(:,:,z)),1) prctile(tovec(image(:,:,z)),99)];
        end

        nexttile([1 1])
        %imshow2(image(:,:,z),cax)
        imagesc(image(:,:,z),cax)
        title(num2str(z))
    end
end
end