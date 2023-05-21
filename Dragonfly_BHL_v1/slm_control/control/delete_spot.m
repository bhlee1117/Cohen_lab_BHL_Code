function delete_spot(app)

idx = cellfun(@str2num,app.SpotListListBox.Value);
app.current_rois(idx)=[];
app.current_xyz(idx,:) = [];
if ~isempty(app.prev_displacement); app.prev_displacement(idx) = []; end
    
pat = wgs_spots(apply_optimal_offsets(app.current_xyz',app.xyz_origin','apply_trans')',30,1);
app.SLM.project(pat)
app.current_mask = pat;