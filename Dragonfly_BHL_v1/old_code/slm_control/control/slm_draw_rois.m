function slm_draw_rois(app)

if isvalid(app.fig_mask_draw)
    figure(app.fig_mask_draw)
else
    app.fig_mask_draw = slm_load_reference_img();
end

hold(app.fig_mask_draw.Children,'on')
cellfun(@(pts) plot(app.fig_mask_draw.Children,pts(:,1),pts(:,2),'o'),app.current_rois)
hold(app.fig_mask_draw.Children,'off')