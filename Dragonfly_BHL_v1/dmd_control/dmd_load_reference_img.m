function fig = dmd_load_reference_img(varargin)
cameraSettings
app = varargin{1};
fig = figure(998);
h_im = findobj(get(fig,'children'),'type','Image');

if isempty(h_im);h_im = imshow2(zeros(500,500),[]);end


if length(varargin) < 2
        fname_im = dir('*pat_ref*'); fname_im = fname_im(end).name;

        % [fname,fpath]=uigetfile('*.*');
        if contains(fname_im,'.fig')
            fig = openfig(fullfile(fpath,fname_im));
        else
            im = imread(fname_im);
            h_im.CData = im;
            app.ref_im = im;
            set(gca,'clim',[min(im(:)) max(im(:))])
            xlim([0 size(im,2)]);ylim([0 size(im,1)]);
        end
        try
            fname = dir('log*'); fname = fname.name;
            fid = fopen(fname);
            snap_list = textscan(fid,'%s (%4d,%4d) %4dx%4d %d %*s %*s %4dx%4d %*s %*s %*s %*s %*s %*s %*s'...
                ,'delimiter',{'\n',' '},'headerlines',1,'multipledelimsasone',true);
    %         snap_list = reshape(snap_list{1},16,[])';
            fclose(fid);
            cam_app = evalin('base','hCameraApp');
            sensor_h = 2048;%cam_app.sensor_size(2);
            sensor_w = 2048;%cam_app.sensor_size(1);
            for ii=1:size(snap_list{1},1)
                if strcmp(snap_list{1}{ii},fname_im(1:6))
    %                 AA_str = snap_list{ii,5}; AA_str(strfind(AA_str,'x'))=' ';
                    AA = [snap_list{4}(ii) snap_list{5}(ii)];

                    if snap_list{6}(ii)
%                         app.HorizontalOffsetEditField.Value = double(sensor_height/2-AA(2)/2);
%                         app.VerticalOffsetEditField.Value = double(sensor_width/2-AA(1)/2);
                        app.HorizontalOffsetEditField.Value = double(sensor_height/2-AA(1)/2);
                        app.VerticalOffsetEditField.Value = double(sensor_width/2-AA(2)/2);
                    else
                        app.HorizontalOffsetEditField.Value = double(snap_list{3}(ii));
                        app.VerticalOffsetEditField.Value = double(snap_list{2}(ii));
                    end
                end
            end
        end
else
        fname_im = dir('*dmd_reg.tif');
        [~,fidx] = max([fname_im.datenum]);
        fname_im = fname_im(fidx).name;
%         pause(0.05)
%         [fname_im,~]=uigetfile('*.*');
%         drawnow; pause(0.05);
        if contains(fname_im,'.fig')
            openfig(fullfile(fname_im));
        else
            im = imread(fname_im);
            imshow(im,[],'initialmagnification','fit')
        end
end


if ~isvalid(app.h_current_rois)
    app.h_current_rois = line(nan,nan,'linewidth',1); 
end
if ~isempty(app.current_rois)
    app.h_current_rois.XData = app.current_rois{:}(:,1);
    app.h_current_rois.YData = app.current_rois{:}(:,2);
    app.h_current_rois.Color = 'r';
else
    app.h_current_rois.XData = nan; 
    app.h_current_rois.YData = nan; 
end

