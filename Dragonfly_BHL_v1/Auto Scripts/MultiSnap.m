n_snaps = 6;

w_time = 10; %min

for ii=1:n_snaps
    hDragonflyApp.SnapTextArea.Value = sprintf('spot_drift_%g',ii);
    str_now = datestr(now,'HHMMSS');
    hDaq.orangeShutter = 1;
    mov = hDragonflyApp.snapOnePicture;    
    imwrite(mov,fullfile([str_now hDragonflyApp.SnapTextArea.Value{:} '.tif']))
    hDaq.orangeShutter = 0;
    
    if ii~=n_snaps
    pause(w_time*60);
    end
end
%%
im_all = [];
[fname,~]=uigetfile('*.tif','multiselect','on');
for ii=1:n_snaps
    im_all(:,:,ii) = imread(fname{ii});
end