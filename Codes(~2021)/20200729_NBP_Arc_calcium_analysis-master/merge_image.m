function merge_image(im,colors,ticks)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[50 50 900 1200]);
[d1 d2 d3]=cellfun(@size,im);
if size(unique(d1),2)>1 || size(unique(d2),2)>1
RGB=zeros(max(d1),max(d2),3);    
else
RGB=zeros(size(im{1},1),size(im{1},2),3);
end

for i=1:size(im,2)
    max_im{i}=double(max(im{i},[],3));
    max_im{i}=max_im{i}./max(max(max_im{i}));
    for j=1:3
        RGB(1:d1(1,i),1:d2(1,i),j)=RGB(1:d1(1,i),1:d2(1,i),j)+max_im{i}.*colors(i,j);
    end
end
imagesc(RGB)
axis equal tight off
if ~isempty(ticks)
    for i=1:size(im,2)
        text(10,25*i,ticks{1,i},'color',colors(i,:),'FontName','arial rounded mt bold','FontSize',10)
    end
end
end