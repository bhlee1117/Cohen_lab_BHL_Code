function slm_draw_rois(app)
df_app = evalin('base','hDragonflyApp');

% if isvalid(app.fig_mask_draw)
%     figure(app.fig_mask_draw)
% end
% figure(app.fig_mask_draw);
app.fig_mask_draw;
imshow(df_app.camApp.snapOneFrame,[],'initialmagnification','fit')

hold(app.fig_mask_draw.Children,'on')
if length(app.current_rois{1}(:,1))==1
    cellfun(@(pts) plot(app.fig_mask_draw.Children,pts(:,1),pts(:,2),'o','linewidth',2),app.current_rois)
else
    cellfun(@(pts) plot(app.fig_mask_draw.Children,pts(:,1),pts(:,2)),app.current_rois)
end
    
hold(app.fig_mask_draw.Children,'off')