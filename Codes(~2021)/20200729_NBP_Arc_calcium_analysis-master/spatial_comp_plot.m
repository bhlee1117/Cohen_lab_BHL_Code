function spatial_comp_plot(A_or,options,cmap,scale)

figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 1200 1200]);

RGB=zeros(options.d1,options.d2,3);
for i=1:size(A_or,2)
    if isempty(cmap)
    color=rand(1,3);
    while sum(color)<1
     color=rand(1,3);   
    end
    else
        color=cmap;
    end
    I=reshape(A_or(:,i),options.d1,options.d2);
    for j=1:3
        RGB(:,:,j)=RGB(:,:,j)+color(1,j)*I*10;
    end
end
imagesc(imresize(RGB,scale))
hold all
axis equal tight
end
