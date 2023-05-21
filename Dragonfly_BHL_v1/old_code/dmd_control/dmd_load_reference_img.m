function fig = dmd_load_reference_img(varargin)
fig = figure(998);clf;
if length(varargin) < 2
        fname_im = dir('*pat_ref*'); fname_im = fname_im(end).name;

        % [fname,fpath]=uigetfile('*.*');
        if contains(fname_im,'.fig')
            fig = openfig(fullfile(fpath,fname_im));
        else
            im = imread(fname_im);
            imshow(im,[],'initialmagnification','fit')
        end
        
        fname = dir('log*'); fname = fname.name;
        fid = fopen(fname);
        snap_list = textscan(fid,'%s (%4d,%4d) %4dx%4d %d %*s %*s %4dx%4d %*s %*s %*s %*s %*s %*s %*s'...
            ,'delimiter',{'\n',' '},'headerlines',1,'multipledelimsasone',true);
%         snap_list = reshape(snap_list{1},16,[])';
        fclose(fid);
        for ii=1:size(snap_list{1},1)
            if strcmp(snap_list{1}{ii},fname_im(1:6))
%                 AA_str = snap_list{ii,5}; AA_str(strfind(AA_str,'x'))=' ';
                AA = [snap_list{4}(ii) snap_list{5}(ii)];

                if snap_list{6}(ii)
                    varargin{1}.HorizontalOffsetEditField.Value = double(1024-AA(1)/2);
                    varargin{1}.VerticalOffsetEditField.Value = double(1024-AA(2)/2);
                else
                    varargin{1}.HorizontalOffsetEditField.Value = double(snap_list{2}(ii));
                    varargin{1}.VerticalOffsetEditField.Value = double(snap_list{3}(ii));
                end
            end
        end
else
%         fname_im = dir('*dmd_ref*').name;
        [fname_im,~]=uigetfile('*.*');
        if contains(fname_im,'.fig')
            openfig(fullfile(fname_im));
        else
            im = imread(fname_im);
            imshow(im,[],'initialmagnification','fit')
        end
end

