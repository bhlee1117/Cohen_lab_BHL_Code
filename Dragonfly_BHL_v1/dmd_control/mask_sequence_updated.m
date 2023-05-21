function mask_sequence_updated(~,evt)

app = evt.AffectedObject;

try
    h_roi_data = cell2mat(cellfun(@(x) [x;nan nan],app.mask_sequence_rois','uniformoutput',false));

    figure(app.fig_mask_draw)
    if ~isvalid(app.h_mask_sequence_rois)
        app.h_mask_sequence_rois = line(nan,nan,'linewidth',1); 
    end
    if ~isempty(h_roi_data)
        app.h_mask_sequence_rois.XData = h_roi_data(:,1);
        app.h_mask_sequence_rois.YData = h_roi_data(:,2);
        app.h_mask_sequence_rois.Color = 'g';
    else
        app.h_mask_sequence_rois.XData = nan; 
        app.h_mask_sequence_rois.YData = nan; 
    end


    roi_label = cellfun(@num2str,num2cell(1:size(app.mask_sequence,3),1),'uniformoutput',false);
    roi_label_pos = cellfun(@(x) mean(x,1,'omitnan')+[10 10],app.mask_sequence_rois,'uniformoutput',false);
    if ~isvalid(app.h_mask_sequence_text)
        app.h_mask_sequence_text = text(zeros(1,length(roi_label)),zeros(1,length(roi_label)),'');
    elseif length(app.h_mask_sequence_text)~=length(roi_label)
        delete(app.h_mask_sequence_text)
        app.h_mask_sequence_text = text(zeros(1,length(roi_label)),zeros(1,length(roi_label)),'');
    end

    if ~isempty(app.h_mask_sequence_text)
        arrayfun(@(h,label,pos) set(h,'String',label{:},'Position',pos{:},'color','g','fontsize',16),...
        app.h_mask_sequence_text ,roi_label',roi_label_pos')
    end
catch
    app.h_mask_sequence_rois = line(nan,nan,'linewidth',1);
    delete(app.h_mask_sequence_text)
    app.mask_sequence_rois = [];
    warning('mask sequence can not be converted to ROI')
end
app.MaskNoTotalEditField.Value = size(app.mask_sequence,3);