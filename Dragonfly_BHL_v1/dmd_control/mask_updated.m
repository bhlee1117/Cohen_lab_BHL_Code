function mask_updated(~,evt)

app = evt.AffectedObject;

h_roi_data = cell2mat(cellfun(@(x) [x;nan nan],app.current_rois','uniformoutput',false));

figure(app.fig_mask_draw)
if ~isvalid(app.h_current_rois)
    app.h_current_rois = line(nan,nan,'linewidth',1); 
end
if ~isempty(h_roi_data)
    app.h_current_rois.XData = h_roi_data(:,1);
    app.h_current_rois.YData = h_roi_data(:,2);
    app.h_current_rois.Color = 'r';
else
    app.h_current_rois.XData = nan; 
    app.h_current_rois.YData = nan; 
end