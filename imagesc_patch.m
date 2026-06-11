function imagesc_patch(image,cax)


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
        im2show=image(:,:,z);
        if nargin<2
            cax=[prctile(tovec(im2show),1) prctile(tovec(im2show),99)];
            if range(cax)==0
                cax=[min(im2show(:)) max(im2show(:))];
            end
        end

        nexttile([1 1])
        %imshow2(image(:,:,z),cax)
        imagesc(image(:,:,z),cax)
        title(num2str(z))
    end
end
end