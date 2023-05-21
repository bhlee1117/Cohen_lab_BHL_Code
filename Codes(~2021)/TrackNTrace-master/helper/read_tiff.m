function imgdata = read_tiff(filename, convert2double, minmax_frame)
% Usage: imgdata = read_tiff(filename, convert2double, minmax_frame)
%
% Loads an image or image stack from a Tiff file. It possible to read only
% a specified part of the file. File can be read in their original format
% (instead of converting to double).
%
% Note: For tif files larger than 4 GB (which are not allowed in the TIFF
% standard but for example written by imageJ) the function automatically 
% switches to read_BigTiff.m
%
% Input:
%   filename: Tiff file to read
%   convert2double: Datatype conversion to double after reading | default: true
%   minmax_frame: [minframe,maxframe] to read. If only one number is specified,
%                 it is assumed to be maxframe and [1, maxframe] is read.
%                 If maxframe==inf the file is read from minframe to the end.| default: read whole movie
%
% NOTE: Deactivate double conversion if not neccessary, reading the same datatype used in the Tiff file.
%       This potentially saves a lot of memory (factor 4 for 16-bit Tiff movies) and thus increased performance.
%       You can still perform a cast to double (or other datatypes) later if desired.
%         Example: imgdata = double( read_tiff(filename) );

% By
% Simon Christoph Stein - August 2015
% E-Mail: scstein@phys.uni-goettingen.de
%

fprintf('Reading file ''%s''.. \n', filename)

if nargin < 2 || isempty(convert2double)
    convert2double = true;
end

if nargin < 3
    minmax_frame = [];
end

% -- Create Tiff object --
% Check filename, append '.tif' if neccessary
[pathstr,name,ext] = fileparts(filename);
if strcmp(ext, '.tif') == 0
    filename = [filename '.tif'];
end

% Check if file is larger then 3.9GB
fileinfo = dir(filename);
if isempty(fileinfo)
   error('Could not open file ''%s''', filename);
end
if( (fileinfo.bytes/(1024*1024*1024) ) > 3.9)
    warning off backtrace
    warning('Reading image file > 3.9GB. Switching to read_BigTiff (unstable).');
    warning on backtrace
    imgdata = read_BigTiff(filename,convert2double,minmax_frame);
    return;
end

try
    warning('off','all');
    t = Tiff(filename,'r');
catch err
    warning('on','all');
    rethrow(err);
end
cleanupTrigger = onCleanup(@() cleanupFunc(t));


% -- Extract header information --
rows = t.getTag('ImageLength'); % 1
cols = t.getTag('ImageWidth'); % 2
% t.getTag('ImageDepth'); % what is this for?

bitdepth = t.getTag('BitsPerSample');
sampleFormat = t.getTag('SampleFormat');
samplesPerPixel = t.getTag('SamplesPerPixel'); % TODO: For multichannel Tiff. Reading not supported here atm.
if samplesPerPixel > 1
    warning(' Multichannel Tiff not supported! Reading first channel..');
end


% -- Count number of frames for pre-allocation --
warning('on','all');
not(t.lastDirectory); % this forces the Tiff object to read the first frames tag data and check for errors.
warning('off','all');

% If no interval was given, determine number of frames and read whole file
if isempty(minmax_frame)
        % Count number of frames in tiff file
        while not(t.lastDirectory)
            t.nextDirectory;
        end
        last_frame = t.currentDirectory;
        minmax_frame = [1, last_frame];
% If one number was given, treat as maxframe
elseif numel(minmax_frame)==1 
        minmax_frame = [1,minmax_frame];
% If maxframe is 0, start at minframe and count to end of file
elseif (numel(minmax_frame)==2) && (minmax_frame(2) == inf)
    try
        t.setDirectory(minmax_frame(1));
    catch err
        error('Minimum frame cannot be read. Movie shorter than minframe?')
    end
    
    % Count number of frames from minframe to end of file
    while not(t.lastDirectory)
        t.nextDirectory;
    end
    last_frame = t.currentDirectory;
    minmax_frame = [minmax_frame(1), last_frame];
end
% nr_frames
nr_frames = minmax_frame(2)-minmax_frame(1)+1;


% -- pre allocate memory according to datatype --
switch sampleFormat
    case Tiff.SampleFormat.IEEEFP
        imgdata = single(zeros(rows, cols, nr_frames));
        sampleFormatName = 'float';
    case Tiff.SampleFormat.UInt
        switch bitdepth
            case 8
                imgdata = uint8(zeros(rows, cols, nr_frames));
            case 16
                imgdata = uint16(zeros(rows, cols, nr_frames));
            case 32
                imgdata = uint32(zeros(rows, cols, nr_frames));
            case 64
                imgdata = uint64(zeros(rows, cols, nr_frames));
            otherwise
                error('Unsupported bitdepth. Supported are 8/16/32/64')
        end
        sampleFormatName = 'uint';
    case Tiff.SampleFormat.Int
        switch bitdepth
            case 8
                imgdata = int8(zeros(rows, cols, nr_frames));
            case 16
                imgdata = int16(zeros(rows, cols, nr_frames));
            case 32
                imgdata = int32(zeros(rows, cols, nr_frames));
            case 64
                imgdata = int64(zeros(rows, cols, nr_frames));
            otherwise
                error('Unsupported bitdepth. Supported are 8/16/32/64')
        end
        sampleFormatName = 'int';
    otherwise
        error('Unsupported datatype ''%s\''. Supported: ''float'', ''uint'', ''int'' ', datatype);
end


% -- File information --
spacing = '    ';
fprintf('%sRows: %i\n',spacing,rows)
fprintf('%sColumns: %i\n',spacing,cols)
fprintf('%sMax Frame: %i\n',spacing,nr_frames)
fprintf('%sBitdepth: %i\n',spacing,bitdepth)
fprintf('%sDatatype: %s\n',spacing,sampleFormatName)


%  -- Rewind to the first frame that should be read --
try
   t.setDirectory(minmax_frame(1));
catch err
    error('Minimum frame cannot be read. Movie shorter than minframe?')
end
msgAccumulator = ''; % needed for one-line printing

% -- read rest of the stack --
startTime = tic;
lastElapsedTime = 0;

for iFrame = 1:nr_frames
    elapsedTime = toc(startTime);
    if( (elapsedTime-lastElapsedTime)>0.25) % output every 0.25 seconds
        rewindMessages()
        rewPrintf('\nReading Frame %i/%i ..',iFrame, nr_frames);
        lastElapsedTime = elapsedTime;
    end
    imgdata(:,:,iFrame) = t.read();
    
    % Check if the last frame has been reached
    if (t.lastDirectory && iFrame ~= nr_frames)
       nr_frames = iFrame;
       imgdata = imgdata(:,:,1:nr_frames); 
       
       % warning output
       warning('on','all');       
       fprintf('\n >> WARNING <<: File has less frames than specified. Read %i frames.\n',nr_frames);
       rewPrintf('\n >> WARNING <<: File has less frames than specified. Read %i frames.\n',nr_frames); % output twice to prevent warning deleted by rewindMessages()               
       break;
    end
    
    if ~(iFrame == nr_frames)
        t.nextDirectory;
    end
end
rewindMessages()
rewPrintf('\nReading Frame %i/%i .. done\n',nr_frames, nr_frames);

if convert2double
    imgdata = double(imgdata);
end

    function rewPrintf(msg, varargin)
        % Rewindable message printing: Print msg and cache it.
        % Usage is analogous to sprintf.
        msg = sprintf(msg, varargin{:});
        msgAccumulator = [msgAccumulator, msg];
        fprintf(msg);
    end

    function rewindMessages()
        % Remove cached messages from command line, reset cache
        reverseStr = repmat(sprintf('\b'), 1, length(msgAccumulator));
        fprintf([reverseStr]);
        
        msgAccumulator = '';
    end
end



function cleanupFunc(t)
%     fprintf('\n')
t.close;
warning('on','all'); % Restore warnings
%     fprintf('\nUSER EXIT');
end