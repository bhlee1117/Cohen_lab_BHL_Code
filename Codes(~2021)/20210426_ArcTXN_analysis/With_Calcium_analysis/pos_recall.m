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