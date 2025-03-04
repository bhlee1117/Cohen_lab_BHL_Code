function current_xyz_updated(~,evt)
app = evt.AffectedObject;
app.SpotListListBox.Items = cellfun(@num2str,num2cell(1:size(app.current_xyz,1),1),'uniformoutput',false);

h_roi_data = cell2mat(cellfun(@(x) [x;nan nan],app.current_rois,'uniformoutput',false));

figure(app.fig_mask_draw)
if ~isvalid(app.h_rois)
    app.h_rois = line(nan,nan,'linewidth',1); 
end
if ~isempty(h_roi_data)
    app.h_rois.XData = h_roi_data(:,1);
    app.h_rois.YData = h_roi_data(:,2);
    app.h_rois.Marker= '*';
    app.h_rois.Color = 'r';
end

roi_label = cellfun(@num2str,app.SpotListListBox.Items','uniformoutput',false);
roi_label_pos = cellfun(@(x) mean(x,1)+[3 3],app.current_rois,'uniformoutput',false);
if ~isvalid(app.h_text)
    app.h_text = text(zeros(1,length(roi_label)),zeros(1,length(roi_label)),'');
elseif length(app.h_text)~=length(roi_label)
    delete(app.h_text)
    app.h_text = text(zeros(1,length(roi_label)),zeros(1,length(roi_label)),'');
end

if ~isempty(app.h_text)
    arrayfun(@(h,label,pos) set(h,'String',label{:},'Position',pos{:},'color','r','fontsize',16),...
    app.h_text,roi_label,roi_label_pos)
end
% app.h_text = text(zeros(1,length(roi_label)),zeros(1,length(roi_label)),'');



