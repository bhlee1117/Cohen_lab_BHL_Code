    function im_recall(handles,im_TXN)
        slider_value = round(get(handles.slider,'Value'));      
        handles.Image=im_TXN(:,:,slider_value);
        imagesc(handles.Image)
        hold all
       title(num2str(slider_value))
        axis equal
axis off
    end