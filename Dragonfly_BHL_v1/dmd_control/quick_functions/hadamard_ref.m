function hadamard_ref(app,m_q)

if ~exist('m_q','var')
    m_q = [63 14];
end
folder = 'E:\';
cal_file_prefix = ['had_cal_1x'];
ref_file_prefix = ['M-YQ0201-12_FOV1_20um_had'];
app.ref_im = dmd_hadamard_reconstruct(app,m_q,folder,cal_file_prefix,ref_file_prefix);
fig = figure(998);
h_im = findobj(get(fig,'children'),'type','Image');
h_im.CData = app.ref_im;
set(gca,'clim',[prctile(app.ref_im(:),1),prctile(app.ref_im(:),99.95)])
end