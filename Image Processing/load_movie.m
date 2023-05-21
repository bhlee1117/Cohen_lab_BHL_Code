function [ims, varargout] = load_movie(fname,varargin)
%LOAD_MOVIE Loads movie dataset in various formats used at the Cohen Lab
%   load_movie reads from dcimg, binary and tif image formats. The first
%   input argument contains a path to read from, either a filename ending
%   in 'dcimg' or 'bin', or a folder name to read tif files from. It
%   optionally accepts a specified list of frame indexes to read. If
%   omitted or empty, all available frames are returned. For binary files,
%   the number of rows and columns must be provided as inputs. For dcimg
%   and tif, the dimensions are not needed. A vectorized movie matrix is
%   returned, with dimensions [nrows*ncols nframes].
%
%   movie = loadmovie('c:\file.bin',frames,nrows,ncols) reads specified
%   frames from a binary file.
%
%   movie = loadmovie('c:\file.bin',nrows,ncols) reads all frames from a
%   binary file. 
%
%   movie = loadmovie('c:\file.dcimg',frames) reads specified frames from a
%   dcimg file.
%
%   movie = loadmovie('c:\file.dcimg') reads all frames from a dcimg file. 
%
%   movie = loadmovie('c:\folder',frames) reads specified frame indexes
%   from a series of tif files listed in the folder.
%
%   movie = loadmovie('c:\folder') reads all frames from a folder with
%   tifs.
%
%   [movie nrows ncols] = loadmovie(...) also returns frame dimensions,
%   useful when reading from dcimg or tif files.
%
%   This function needs readDCIM to read dcimg files, and readBinMov to
%   read binary files.
%
%   2014 Vicente Parot
%   Cohen Lab - Harvard University
%

%% check input and output format, get data format
switch nargin
    case 0
        error('missing read path')
    case 1
        frames_to_read = [];
    case 2
        frames_to_read = varargin{1};
    case 3
        frames_to_read = [];
        nr = varargin{1};
        nc = varargin{2};
    case 4
        frames_to_read = varargin{1};
        nr = varargin{2};
        nc = varargin{3};
    otherwise
        error('wrong input arguments format')
end
if nargout > 3
    error('wrong output arguments format')
end
if strcmp(fname(max(1,end-4):end),'dcimg')
    % dcimg -> ims
    if nargin > 2
        warning('ignoring image dimension input arguments')
    end
    if nargout < 3
        warning('would you like to have the image dimensions? call with three output arguments')
    end
    data_format = 'dcimg';
elseif strcmp(fname(max(1,end-2):end),'bin')
    % bin -> ims
    if nargin < 3
        error('please provide image dimensions')
    elseif nargout > 1
        warning('passing dimension inputs as dimension outputs')
    end
    data_format = 'bin';
else
    % \*.tif -> ims
    if nargin > 2
        warning('ignoring image dimension input arguments')
    end
    if nargout < 3
        warning('would you like to have the image dimensions? call with three output arguments')
    end
    data_format = 'tif';
end

%% load files
% figure out dcimg or tif folder
switch data_format
    case 'dcimg'
        % dcimg -> ims
        evalc('[~,nframes] = dcimgmatlab(1, fname);');
        if ~nframes %#ok nframes is defined in evalc call
            error('no frames in dcimg file')
        end
        nframes = double(nframes);
        if isempty(frames_to_read)
            frames_to_read = 1:nframes; %#ok frames_to_read will be used in an evalc call
        end
        tic
        disp 'reading dcimg file ...'
        evalc('ims = readDCIMG(fname,frames_to_read);'); % supress readDCIMG output
        disp(['readDCIMG took ' num2str(toc) ' s']);
    case 'bin'
        % bin -> ims
        tic
        disp 'reading binary file ...'
%         fileProps = dir(fname);
%         nFrames = fileProps.bytes/(nr*nc*2);
        if ~isempty(frames_to_read)
            ims = readBinMov_times(fname,nr,nc, frames_to_read);
        else
            ims = readBinMov(fname, nr, nc);
        end;
        disp(['reading took ' num2str(toc) ' s']);
    case 'tif'
        % \*.tif -> ims
        fname = [fname filesep]; % append a filesep, no harm
        filenames = dir([fname '*tif']);
        nframes = numel(filenames);
        lastimage = imread([fname filenames(end).name]);
        [nr,nc] = size(lastimage);
        if isempty(frames_to_read)
            frames_to_read = 1:nframes;
        end
        nt = numel(frames_to_read);
        ims = zeros(nr,nc,nt,class(lastimage));
        tic
        disp 'reading tif files ...'
        for it = 1:nt
            ims(:,:,it) = imread([fname filenames(frames_to_read(it)).name]);
            if ~mod(it,250), disp([num2str(round(it/nt*100)) '% read in ' num2str(toc) ' s']); end
        end
        disp([num2str(nt) ' files read in ' num2str(toc) ' s']);
    otherwise
        why
end
[nr, nc, ~] = size(ims);
tovec = @(m)reshape(m,nr*nc,[]);
ims = tovec(ims);
ims = double(ims);
varargout{1} = nr;
varargout{2} = nc;
