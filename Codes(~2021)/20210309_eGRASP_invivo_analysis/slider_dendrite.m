%function slider_dendrite(candidatePos,D_coor,c_image_xy)
for dend=1:size(D_coor{1,1},2)
rng=round([min(D_coor{1,1}{dend}); max(D_coor{1,1}{dend})]);
[crop_im{1}]=c_image_xy{1}(rng(1,2):rng(2,2),rng(1,1):rng(2,1),:);
s=[candidatePos{1,1}(dend,3)-candidatePos{1,1}(dend,1);
     candidatePos{1,1}(dend,4)-candidatePos{1,1}(dend,2)];
[crop_im{2}]=c_image_xy{2}(rng(1,2)+s(2):rng(2,2)+s(2),rng(1,1)+s(1):rng(2,1)+s(1),:);
end
handles.fig = figure;
handles.axes1 = axes('Units','pixels','Position',[50 50 400 400]);
handles.Image = crop_im;

handles.fig = figure;
handles.axes2 = axes('Units','pixels','Position',[50 50 400 400]);
handles.Image = imCal;


handles.cell=cell_plane;
handles.resize_ratio=1;
%rsz=size(max_res_im,1)/size(im_ca_crop,1);
imagesc(handles.Image);
colormap('gray')
axis equal
axis off
hold all

plot(cell_plane(:,1)/1,cell_plane(:,2)/1,'ro','markersize',13)

handles.slider = uicontrol('Style','slider','Position',[55 20 500 20],'Min',1,'Max',size(candidatePos,1),'Value',1);
handles.Listener = addlistener(handles.slider,'Value','PostSet',@(s,e) pos_recall(handles,cell_list,zmat));
guidata(handles.fig);
end

    function pos_recall(handles,cell_list,zmat)
        slider_value = round(get(handles.slider,'Value'));      
       cell_plane=cell_list(find(round(cell_list(:,3))>=zmat-slider_value & round(cell_list(:,3))<=zmat+slider_value),1:3);
       resize_ratio=handles.resize_ratio;
       imagesc(handles.Image) 
        axis equal
        colormap('gray')
       hold all
       plot(cell_plane(:,1)/resize_ratio,cell_plane(:,2)/resize_ratio,'ro','markersize',13)
       text(3,10,['Zmat=',num2str(zmat),'   ',num2str(zmat-slider_value),'~',num2str(zmat+slider_value),'    Slider=',num2str(slider_value)],'color',[1 1 1])
       
axis off
    end