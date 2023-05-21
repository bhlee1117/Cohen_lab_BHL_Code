function dmd_draw_rois(app)
df_app = evalin('base','hDragonflyApp');

if isvalid(app.fig_mask_draw)
    figure(app.fig_mask_draw)
end

app.fig_mask_draw;
im = findobj(app.fig_mask_draw.Children,'type','image');
im.CData = df_app.camApp.snapOneFrame;
set(gca,'clim',[min(im.CData(:)) max(im.CData(:))])
% imshow(df_app.camApp.snapOneFrame,[],'initialmagnification','fit')

% hold(app.fig_mask_draw.Children,'on')
% cellfun(@(pts) plot(app.fig_mask_draw.Children,pts(:,1),pts(:,2)),app.current_rois)
% hold(app.fig_mask_draw.Children,'off')