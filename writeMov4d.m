function writeMov4d(filename, movie, opt)
% writeMov4d(filename, movie, opt)
%
% Writes a 4-D movie (X x Y x Ch x T) to an MPEG-4 file.
%
% REQUIRED
%   filename  : output path (string), e.g. 'myMovie.mp4'
%   movie     : X x Y x Ch x T numeric array
%
% opt FIELDS (all optional — defaults shown)
%   opt.frameRate      : playback frame rate of the output video  [10]
%   opt.timeVec        : 1 x T vector; frame timestamp in ms shown as title
%                        If empty, frame index is shown instead.          []
%   opt.range          : [lo hi] display clim. Auto from mid-frame if []  []
%   opt.roiXY          : N x 2 [x y] coordinates of ROI overlay          []
%   opt.labelVec       : 1 x T logical / binary vector; frames where
%                        labelVec==1 get a text label in the top-left      []
%   opt.labelText      : string shown when labelVec==1              ['stim']
%   opt.labelColor     : text color for label                       [1 1 0]
%   opt.colormap       : colormap name or M x 3 matrix              'turbo'
%   opt.colorbarLabel  : colorbar label string                    'Z score'
%   opt.quality        : VideoWriter quality (1-100)                  [100]
%   opt.scaleBar       : length of scale bar in pixels. No bar if []      []
%   opt.scaleBarText   : label drawn above the bar (e.g. '50 µm')        ''
%   opt.scaleBarColor  : color of bar and text                      [1 1 1]
%   opt.scaleBarWidth  : line width of the bar                           [3]
%   opt.fontSize       : font size applied to all axes / labels           []
%   opt.fontName       : font name applied to all axes / labels           []
%
% EXAMPLE
%   opt.frameRate    = 30;
%   opt.timeVec      = linspace(0, 500, size(movie,4));
%   opt.range        = [-3 3];
%   opt.labelVec     = stimOn;   % logical 1 x T
%   opt.labelText    = 'stim';
%   writeMov4d('out.mp4', movie, opt);

% ---------- defaults ----------
if nargin < 3, opt = struct(); end
opt = setdefault(opt, 'frameRate',    10);
opt = setdefault(opt, 'timeVec',      []);
opt = setdefault(opt, 'range',        []);
opt = setdefault(opt, 'roiXY',        []);
opt = setdefault(opt, 'labelVec',     []);
opt = setdefault(opt, 'labelText',    'stim');
opt = setdefault(opt, 'labelColor',   [1 1 0]);
opt = setdefault(opt, 'colormap',     'turbo');
opt = setdefault(opt, 'colorbarLabel','Z score');
opt = setdefault(opt, 'quality',      100);
opt = setdefault(opt, 'scaleBar',     []);
opt = setdefault(opt, 'scaleBarText', '');
opt = setdefault(opt, 'scaleBarColor',[1 1 1]);
opt = setdefault(opt, 'scaleBarWidth',3);
opt = setdefault(opt, 'fontSize',     []);
opt = setdefault(opt, 'fontName',     []);

nFrames = size(movie, 4);

% auto display range from middle frame
if isempty(opt.range)
    midSlice = squeeze(movie(:,:,:, floor(nFrames/2)));
    opt.range = [min(midSlice(:)) max(midSlice(:))];
end

% validate labelVec
useLabelVec = ~isempty(opt.labelVec);
if useLabelVec
    opt.labelVec = logical(opt.labelVec(:)');   % ensure 1 x T logical
    if numel(opt.labelVec) ~= nFrames
        warning('writeMov4d: labelVec length (%d) != nFrames (%d). Ignoring.', ...
            numel(opt.labelVec), nFrames);
        useLabelVec = false;
    end
end

% ---------- open video ----------
set(gcf, 'Color', 'k');

myVideo = VideoWriter(filename, 'MPEG-4');
myVideo.FrameRate = opt.frameRate;
myVideo.Quality   = opt.quality;
open(myVideo);

% ---------- write frames ----------
for i = 1:nFrames
    clf;
    pbaspect([size(movie,2) size(movie,1) 1] * 2);

    h = imshow2(squeeze(movie(:,:,:,i)), opt.range);   %#ok<NASGU>
    axis tight off equal;
    hold on;

    % ROI overlay
    if ~isempty(opt.roiXY)
        plot(opt.roiXY(:,1), opt.roiXY(:,2), ...
            'Color', [0 0.6 1], 'LineWidth', 2);
    end

    % title: timestamp or frame index
    if ~isempty(opt.timeVec)
        title([num2str(opt.timeVec(i)) ' ms'], 'Color', [1 1 1]);
    else
        title(['Frame ' num2str(i)], 'Color', [1 1 1]);
    end

    % top-left label for flagged frames
    if useLabelVec && opt.labelVec(i)
        text(0.02, 0.97, opt.labelText, ...
            'Units', 'normalized', ...
            'Color', opt.labelColor, ...
            'FontSize', 14, ...
            'FontWeight', 'bold', ...
            'VerticalAlignment', 'top');
    end

    % scale bar
    if ~isempty(opt.scaleBar)
        hBar = drawScaleBar(opt.scaleBar, 'horizontal', ...
            'Color',     opt.scaleBarColor, ...
            'LineWidth', opt.scaleBarWidth);
        if ~isempty(opt.scaleBarText)
            % place text centred above the bar
            xData = get(hBar, 'XData');
            yData = get(hBar, 'YData');
            text(mean(xData), yData(1)-2, opt.scaleBarText, ...
                'Color',               opt.scaleBarColor, ...
                'FontSize',            11, ...
                'FontWeight',          'bold', ...
                'HorizontalAlignment', 'center', ...
                'Fontname',opt.fontName,...
                'VerticalAlignment',   'bottom');
        end
    end

    % colormap + colorbar
    colormap(opt.colormap);
    cb = colorbar;
    cb.Label.String = opt.colorbarLabel;
    cb.FontSize  = 15;
    cb.Color     = [1 1 1];

    % font size / name applied to all axes in the figure
    if ~isempty(opt.fontSize), set_fontsize(opt.fontSize); end
    if ~isempty(opt.fontName), set_font(opt.fontName);     end

    frame = getframe(gcf);
    writeVideo(myVideo, frame);
end

close(myVideo);
close(figure(1));
end

% -----------------------------------------------------------------------
function opt = setdefault(opt, field, val)
    if ~isfield(opt, field) || isempty(opt.(field))
        opt.(field) = val;
    end
end
