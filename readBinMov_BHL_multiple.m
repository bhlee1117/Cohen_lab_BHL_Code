function mov = readBinMov_BHL_multiple(filepath, n, time_read, time_segment, prefix)
% readBinMov_BHL_multiple  Read frames from segmented motion-corrected bin files.
%
% Usage:
%   mov = readBinMov_BHL_multiple(filepath)
%   mov = readBinMov_BHL_multiple(filepath, n)
%   mov = readBinMov_BHL_multiple(filepath, n, time_read)
%   mov = readBinMov_BHL_multiple(filepath, n, time_read, time_segment)
%   mov = readBinMov_BHL_multiple(filepath, n, time_read, time_segment, prefix)
%
% Inputs:
%   filepath     - folder containing the bin files and output_data.mat
%   n            - Device_Data index for ROI/size  (default: 3)
%   time_read    - global frame indices to read    (default: all frames)
%   time_segment - frames per bin file             (default: 15000)
%   prefix       - filename prefix before the zero-padded segment number
%                  (default: 'mc_')
%                  Examples:
%                    'mc_'    -> mc_01.bin, mc_02.bin, ...
%                    'mc'     -> mc01.bin,  mc02.bin,  ...
%                    'mc1_0'  -> mc1_001.bin, mc1_002.bin, ...  (auto-detected width)
%                    'mc2'    -> mc201.bin, mc202.bin, ...
%
% Output:
%   mov - H x W x T double array of requested frames

%% ---- defaults ----------------------------------------------------------
if nargin < 2 || isempty(n),            n            = 3;     end
if nargin < 4 || isempty(time_segment), time_segment = 15000; end
if nargin < 5 || isempty(prefix),       prefix       = 'mc_'; end

%% ---- load geometry -----------------------------------------------------
load(fullfile(filepath, 'output_data.mat'), 'Device_Data');
sz = double(Device_Data{n}.ROI([2 4]));   % [nCol nRow]  ->  sz(1)=nCol, sz(2)=nRow

%% ---- default time range ------------------------------------------------
if nargin < 3 || isempty(time_read)
    time_read = 1 : max(Device_Data{1,2}.Counter_Inputs.data);
end
time_read(time_read <= 0) = [];

%% ---- discover matching bin files ---------------------------------------
pattern     = [prefix '*.bin'];
isMCtrace   = dir(fullfile(filepath, pattern));

%% ---- compute which file and local frame each requested frame lives in --
file_read  = ceil(time_read / time_segment);
frame_read = mod(time_read,  time_segment);
frame_read(frame_read == 0) = time_segment;

%% ---- read --------------------------------------------------------------
mov = [];

if ~isempty(isMCtrace)
    % Sort files by name so segment order is correct
    [~, sortIdx] = sort({isMCtrace.name});
    isMCtrace    = isMCtrace(sortIdx);
    file2read= unique(file_read, 'stable');
    if size(file2read,1)==1
    else
        file2read=file2read';
    end

    for j = file2read
        if j > numel(isMCtrace)
            warning('readBinMov_BHL_multiple: segment %d requested but only %d files found matching "%s".', ...
                j, numel(isMCtrace), pattern);
            continue;
        end
        fname = fullfile(filepath, isMCtrace(j).name);
        disp(['Loading ' fname]);
        frames_in_seg = frame_read(file_read == j);
        mov = cat(3, mov, double(readBinMov_times(fname, sz(2), sz(1), frames_in_seg)));
    end

else
    % Fallback: no segmented files found, read raw frames1.bin directly
    warning('readBinMov_BHL_multiple: no files matching "%s" found in %s. Falling back to frames1.bin.', ...
        pattern, filepath);
    mov = double(readBinMov_times(fullfile(filepath, 'frames1.bin'), sz(2), sz(1), time_read));
end

end
