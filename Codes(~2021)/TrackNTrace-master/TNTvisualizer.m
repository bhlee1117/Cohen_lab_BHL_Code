% TrackNTrace: A simple and extendable MATLAB framework for single-molecule localization and tracking
%
%     Copyright (C) 2016  Simon Christoph Stein, scstein@phys.uni-goettingen.de
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function [h_main, movie] = TNTvisualizer(movieOrTNTfile, candidateDataOrTNTfile, candidateParams, refinementData, refinementParams, trackingData, trackingParams, FPS, is_blocking, firstFrame)
% Visualizer for movies and TrackNTrace results.
%
% COMMON USAGEs: 
%    (x) TNTvisualizer();
%        Opens a dialog to choose a movie or TNT file from a completed TNT
%        run. When choosing a TNT file, the corresponding movie is loaded
%        automatically along all results from the evaluation.
%    (x) [~,movie] = TNTvisualizer();
%        Same as above, but the loaded movie is returned.
%    (x) [~,movie] = TNTvisualizer(PathToTNTfile);
%        Load TNT file. This loads the corresponding movie and all results from the evaluation.
%    (x) TNTvisualizer(movie);
%        Visualize movie already loaded in MATLAB.
%    (x) [~,movie] = TNTvisualizer(PathToMovie);
%        Same as above, but movie is loaded from disk first.
%    (x) TNTvisualizer(movie, PathToTNTfile);
%        Visualize already loaded movie along with evaluated data from
%        specified TNT file. Make sure the movie fits to the evaluated data.
%    (x) [~,movie] = TNTvisualizer(PathToMovie, PathToTNTfile); 
%        Same as above, but movie is loaded from disk first.

% [ Full USAGE: [GUIhandle, movie] = TNTvisualizer(movieOrTNTfile, candidateDataOrTNTfile, candidateParams, refinementData, refinementParams, trackingData, trackingParams, FPS, is_blocking, firstFrame) ]
%
% Input:
%   movieOrTNTfile: 
%     EITHER: 3D matrix (rows,cols,frames) with intensity values.
%     OR: Filepath to TIF movie -> Loads movie
%     OR: Filepath to TNT file. -> Loads movie and evaluation
%   candidateDataOrTNTfile: 
%     EITHER: 1D cell array, where cell{F} is a 2D matrix with dimensions 
%       N(F)x(2+PC) saving the data of frame F. N(F) is the number
%       of candidates detected in frame F and the columns are 
%       'x','y' plus PC arbitrary additional parameters.
%     OR: Filepath to TNT file. -> Loads evaluated data 
%         (movie must have been given as the first parameter).
%   candidateParams:
%     1D cell array with (2+PC) strings containing the names of the parameters
%     (columns of each cell) of candidateData (above).
%   refinementData: 
%     1D cell array, where cell{F} is a 2D matrix with dimensions
%     N(F)x(3+PF) saving the data of frame F. N(F) is the number
%     of fitted emitters detected in frame F and the columns are 
%     'x','y','z' position plus PF arbitrary additional parameters.
%   refinementParams:
%     1D cell array with (3+PF) strings containing the names of the parameters
%     (columns of each cell) of refinementData (above).
%   trackingData: 
%     2D matrix with dimensions Nx(5+PT). N is the number of all positions
%     tracked positions. The columns are 'TrackID','frame','x','y','z'
%     plus PT arbitrary additional parameters.
%   trackingParams:
%     1D cell array with (5+PT) strings containing the names of the parameters
%     (columns) of trackingData (above).
%   FPS:
%     Frames per second the visualizer plays set to use for playback on
%     startup.
%   is_blocking:
%     If true, the GUI blocks MATLAB execution until it closes.
%   firstFrame:
%     Index of first frame. This is used for displaying the true index
%     (with respect to the full movie) of a frame  if only a certain interval
%     of a movie was processed by TrackNTrace.
%     Example: If frame 201 to 250 from some movie were processed. Setting
%     firstFrame=201 shows 'Fr. 10/50 (210 in movie)' if the 10th processed
%     frame is shown.
%
%  Inputs can be left empty [] for default values.
%
% Output:
%   h_main:
%     Handle to the visualizer figure / GUI.
%   movie:
%     The visualized movie. Useful if movie was loaded by the visualizer.
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2016
%

% Add all paths required to run TNT
addPathsVisualizer();

% Check MATLAB version
MATLABversion = strsplit(version,'.');
if(str2double(MATLABversion(1))>=8 && str2double(MATLABversion(2))>=6) % Use 'drawback nocallbacks' for MATLAB 2015b (8.6.x) and later
    MATLAB_2015b_or_newer = true; 
else
    MATLAB_2015b_or_newer = false; 
end

% Show filechooser dialog to choose movie / TNT file if started without input arguments.
if nargin==0
   [filename, path] = uigetfile({'*.mat;*.tif','Visualizer files'},'Select movie or TNT file to visualize.');
    if isfloat(filename);
        return;
    end; % User clicked cancel
    [~,~,ext] = fileparts(filename);
    switch ext
        case '.tif'
            movieOrTNTfile = [path,filename];
        case '.mat' 
            movieOrTNTfile = [path,filename];
        otherwise
            error('Only movies or TNT files are supported for loading!');
    end
end


% -- Preparing the GUI --
h_main = openfig(['subfun',filesep,'TNTvisualizer_Layout.fig'],'new','invisible');
set(h_main,'Renderer','painters');
set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'CloseRequestFcn',@onAppClose); % Executed on closing for cleanup
set(h_main,'Toolbar','figure');   % Add toolbar needed for zooming
set(h_main, 'DoubleBuffer', 'on') % Helps against flickering
movegui(h_main,'center');

h_all = guihandles(h_main); % Get handles of all GUI objects
axes(h_all.axes); % Select axis for drawing plots

% -- Datacursor --
dcm_obj = datacursormode(h_main);
% defaultDatatipFunction = get(dcm_obj,'UpdateFcn');
set(dcm_obj,'Enable','on');

% -- Setup and input parsing --
% These variables are accessed/modified across different functions of the GUI
allFramesCandidateData = []; % Computed on demand for histogram plots
allFramesRefinementData = []; % Computed on demand for histogram plots

use_bw = false; % Is visualization black/white or colored?
mode = 'movie'; % 'movie' 'candidate'/'refinement'/'tracking'
traj_lifetime = 0; % Trajectories are kept for display for traj_lifetime after the particles last appearance
nr_track_colors = 20; % Size of colorpool to color tracks
traj_displayLength = inf; % Tracking data is displayed for the last traj_displayLength frames only (with respect to the current frame).

movie = []; % Set by parse_inputs_and_setup depending on 'movieOrTNTfile'
candidateData = []; % Set by parse_inputs_and_setup depending on 'candidateDataOrTNTfile'
id_tracks = []; % Unique ids for the tracks. If tracks are numbered continously from 1 to N, this is identical to the track index.
cell_traj = {}; % Each cell saves the data for one track. Is initialized by parse_inputs_and_setup.
n_tracks = 0;
parse_inputs_and_setup(nargin);

tracksVisibleInFrame = {}; % tracksVisibleInFrame{frame} is list of all track numbers visible in frame 'frame'
nr_VisibleTracksInFrame = []; % Number of visible tracks in each frame. 
maxNr_visibleTracksInFrame = 0; % Maximum number of tracks concurrently visible in one frame
compute_tracksVisibleInFrame(); % Sets tracksVisibleInFrame, nr_VisibleTracksInFrame and maxNr_visibleTracksInFrame.
% ---------------------------

% Show the GUI!
% (We start invisible in case a movie is loaded during parse_inputs_and_setup)
set(h_main,'visible','on'); 

%  -- Setup UI elements --
% Buttons
set(h_all.button_movieMode,'Callback', {@callback_changeMode,'movie'});
set(h_all.button_candidateMode,'Callback', {@callback_changeMode,'candidate'});
set(h_all.button_refinementMode,'Callback', {@callback_changeMode,'refinement'});
set(h_all.button_trackingMode,'Callback', {@callback_changeMode,'tracking'});

set(h_all.but_play,'Callback',@playCallback);
set(h_all.but_contrast,'Callback',@contrastCallback);
set(h_all.but_autocontrast,'Callback',@autocontrastCallback);
set(h_all.but_autocontrast,'TooltipString',sprintf('Set contrast automatically.\n The used algorithm can be selected in the popup menu to the right.\n Press  and hold "Shift" key during playback to adjust contrast automatically for each frame.'));
set(h_all.popup_autocontrast, 'TooltipString', sprintf('Algorithm used for autocontrast.\n\b\b Spots: Emphasize highest 25%% intensity values.\n\b\b Min/Max: Spans all values.\n\b\b 98%% range: Cuts the lower and upper 1%% of intensities. '));
set(h_all.but_distribution,'Callback',@distributionCallback);

% Slider
if(size(movie,3)>1) 
    set(h_all.slider,'Value',1, 'Min',1,'Max',size(movie,3),'SliderStep',[1/size(movie,3) 1/size(movie,3)],'Callback', @sliderCallback);
else % For single images we disable slider and play button
    set(h_all.slider,'Enable','off');
    set(h_all.but_play,'Enable','off');
end
hLstn = addlistener(h_all.slider,'ContinuousValueChange',@updateSlider); %#ok<NASGU> % Add event listener for continous update of the shown slider value

% Edit fields
set(h_all.edit_FPS,'String',sprintf('%i',FPS), 'Callback', @fpsCallback);
set(h_all.edit_distributionBins, 'Callback', {@callback_intEdit,1,inf});
setNum(h_all.edit_distributionBins, 50, true);
set(h_all.edit_distributionRange,'Callback',{@callback_FloatEdit,0,100});
setNum(h_all.edit_distributionRange,100);

% Checkbox
set(h_all.cb_bw, 'Value', use_bw, 'Callback',@bwCallback);

% Popupmenu
set(h_all.popup_distribution, 'String', 'No data');

%  -- Candidate UI elements --

%  -- Refinement UI elements --

%  -- Tracking UI elements --
set(h_all.edit_lifetime,'String',sprintf('%i',traj_lifetime), 'Callback', @callback_TrajLifetime);
set(h_all.edit_colors,'String',sprintf('%i',nr_track_colors), 'Callback', @callback_trackColors);
set(h_all.edit_trajDisplayLength,'String',sprintf('%i', traj_displayLength), 'Callback', @callback_dispLength);

% -- Timer -> this controls playing the movie --
h_all.timer = timer(...
    'ExecutionMode', 'fixedDelay', ...    % Run timer repeatedly
    'Period', round(1/FPS*1000)/1000, ... % Initial period is 1 sec. Limited to millisecond precision
    'TimerFcn', @onTimerUpdate, ...
    'StartFcn', @onTimerStart, ...
    'StopFcn',  @onTimerStop); % Specify callback

% Store handles to the plot objects (which is faster)
% Handles are set on first use (mostly in plotFrame)
linehandles = -1*ones(maxNr_visibleTracksInFrame,1); % Note: Size can increase dynamically if more lines are needed.
linehandleNr_to_TrackNr = -1*ones(maxNr_visibleTracksInFrame,1); % Stores which track the linehandle is assigned to
dothandle_fit = -1;
dothandle_cand = -1;
imagehandle = -1;

% Draw the marker color depending on background color
track_colors = [];
marker_color = [];
marker_fill_color = [];
drawColors(nr_track_colors);

% -- Variables for playback --
timePerFrame = round(1/FPS*1000)/1000; % limit to millisecond precision
elapsed_time = 0;
frame = 1;

% Plot first frame to get limits right
% Set x,y,color limits
xl = [0.5,size(movie,2)+0.5];
yl = [0.5,size(movie,1)+0.5];
firstImg = movie(:,:,1);
zl = [min(firstImg(:)), max(firstImg(:))];

plotFrame(frame);

xlim(xl);
ylim(yl);
caxis(zl);

updateTopText();

% -- Change into the right mode (candidate/refinement/tracking) --
callback_changeMode();

% For is_blocking==true we stop scripts/functions
% calling the GUI until the figure is closed
if(is_blocking)
    uiwait(h_main);
    drawnow; % makes figure disappear instantly (otherwise it looks like it is existing until script finishes)
end


% --- Nested Functions ---


% Change visualizer into the chosen mode 'movie' 'candidate','refinement','tracking'
% and display its relevant content.
% If "modus" input is given, the mode is set to "modus". Its
% implemented this way to use one callback for all buttons selecting the modes.
    function callback_changeMode(~,~,modus)
        % Note: nargin>2 is true if callback was invoked by a button
        % (and not from a direct call in this file)
        if nargin>2
            % Do nothing if the current modes button is pressed again
            if strcmp(modus,mode)
                return
            end
            mode = modus;
        end
        DEFAULT_COLOR = [0.941,0.941,0.941]; % Default color of buttons.
        SELECTED_COLOR = [0.65, 0.9, 0]; % Color of selected button.
        
        %Reset button colors
        set(h_all.button_movieMode,'BackgroundColor', DEFAULT_COLOR);
        set(h_all.button_candidateMode,'BackgroundColor', DEFAULT_COLOR);
        set(h_all.button_refinementMode,'BackgroundColor', DEFAULT_COLOR);
        set(h_all.button_trackingMode,'BackgroundColor', DEFAULT_COLOR);
        
        % Set all mode specific panels invisible
        set(h_all.panel_tracking,'Visible','off');
        set(h_all.panel_histogram,'Visible','off');
                
        % Delete all graphics objects, except the movie frame and invalidate all handles       
        resetGraphics();
        
        % Delete active datatip
        % WARNING this also deletes other hggroup objects associated with the figure which are also invisible and draggable.
        delete(findall(h_main,'Type','hggroup','HandleVisibility','off','Draggable','on'));
        
        % Mode specific changes (setting datatip function, highlight buttons etc.)
        switch mode
            case 'movie'
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_movieMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', 'No data');
            case 'candidate'
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_candidateMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', candidateParams);
                set(h_all.panel_histogram,'Visible','on');
            case 'refinement'
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_refinementMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', refinementParams);
                set(h_all.panel_histogram,'Visible','on');
            case 'tracking'
                initializeLinehandles(); %Initializes all needed line handles, if this is uncommented, linehandles are created on the fly
                
                set(dcm_obj,'UpdateFcn',{@modeSpecificDatatipFunction});
                set(h_all.button_trackingMode,'BackgroundColor', SELECTED_COLOR);
                set(h_all.popup_distribution, 'String', trackingParams);
                set(h_all.panel_histogram,'Visible','on');
                set(h_all.panel_tracking,'Visible','on');
            otherwise
                error('Unkown mode ''%s''!', mode);
        end
        % Select first histogram entry (There should always be a first entry!)
        % If not set, this runs into problems if selecting entry X in some
        % mode and there are less then X parameters when the mode is switched.
        set(h_all.popup_distribution, 'Value', 1);        
        
        % Resize GUI
        resizeGUIforMode();
        
        % Replot
        updateFrameDisplay();
    end

    % Delete all graphics objects, except the movie frame and invalidate all handles       
    function resetGraphics()
        allHandles = get(h_all.axes,'Children');
        if imagehandle ~= -1 % Remove image handle from list
            allHandles(allHandles==imagehandle) = [];
        end
        delete(allHandles);
        dothandle_fit = -1;
        dothandle_cand = -1;
        linehandles = -1*ones(maxNr_visibleTracksInFrame,1);
    end

% Update size of GUI, show mode specific panels
    function resizeGUIforMode()
        units = 'characters';
        BOTTOM_SPACING = 0.5;
        
        % Get position of last element in GUI (mode dependent)
        switch mode
            case 'movie'
                set(h_all.panel_player,'Units',units);
                pos = get(h_all.panel_player,'Position');
            case 'candidate'
                set(h_all.panel_histogram,'Units',units);
                pos = get(h_all.panel_histogram,'Position');
            case 'refinement'
                set(h_all.panel_histogram,'Units',units);
                pos = get(h_all.panel_histogram,'Position');
            case 'tracking'
                set(h_all.panel_tracking,'Units',units);
                pos = get(h_all.panel_tracking,'Position');
        end
        
        set(h_main,'Units',units);
        win_pos = get(h_main,'Position');
        win_top_pos = win_pos(2)+win_pos(4); % Save top position
        diff_height = pos(2)-BOTTOM_SPACING;
        
        % To resize the figure properly, we first need to move all objects
        % inside.. (Matlab ..)
        all_uiObjects = get(h_main,'Children');
        
        for iObj = 1:numel(all_uiObjects)
            try %#ok<TRYNC>
                set(all_uiObjects(iObj),'Units',units);
                pos = get(all_uiObjects(iObj),'Position');
                pos(2) = pos(2) - diff_height;
                set(all_uiObjects(iObj),'Position',pos);
            end
        end
        
        win_pos(4) = win_pos(4)-diff_height; % Set new window height
        win_pos(2) = win_top_pos-win_pos(4); % Keeps the top position constant
        set(h_main,'Position', win_pos);
        
        % Reset units back to normalized, so figure resizes "properly" (cough..)
        set(h_main,'Units','normalized');
        for iObj = 1:numel(all_uiObjects)
            try %#ok<TRYNC>
                set(all_uiObjects(iObj),'Units','normalized');
            end
        end
    end

    % Custom function for datacursor which shows data relevant to the
    % current mode when clicking the currently plotted data.
    function txt = modeSpecificDatatipFunction(~,event_obj)
       
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        graphObjHandle = get(event_obj,'Target'); % The target object (line/image) of the cursor
        
        if(isgraphics(graphObjHandle,'image')) % Image is selected
            txt = {['X: ',num2str(pos(1))],...
                   ['Y: ',num2str(pos(2))],...
                   ['Intensity: ', num2str(movie(pos(2),pos(1),frame))]};
        else % Plotted position is selected
            I = get(event_obj, 'DataIndex');
            txt = {};
            switch mode
                case 'movie' % -> Image is selected handled above
                case 'candidate'
                    % Plot all parameters available for selected datapoint in the datacursor window
                    for iPar=1:numel(candidateParams)
                        txt = [txt, {[candidateParams{iPar},': ', num2str(candidateData{frame}(I,iPar))]}]; %#ok<AGROW>
                    end
                case 'refinement'
                    % Plot all parameters available for selected datapoint in the datacursor window
                    for iPar=1:numel(refinementParams)
                        txt = [txt, {[refinementParams{iPar},': ', num2str(refinementData{frame}(I,iPar))]}]; %#ok<AGROW>
                    end
                case 'tracking'
                    handleNr = linehandles==graphObjHandle; % Find lineobject for the selected point
                    TrackNr = linehandleNr_to_TrackNr(handleNr);
                    PointData = cell_traj{TrackNr}(I,:); % Data of the selected point
                    TrackID = sprintf('%i',id_tracks(TrackNr)); % Get track ID from its index (in case TracIDs go from 1 to N without missing numbers, TrackNr==TrackID)
                    % Plot all parameters available for selected datapoint in the datacursor window
                    txt = [txt, {['TrackID: ', TrackID]}];
                    for iPar=2:numel(trackingParams)
                        txt = [txt, {[trackingParams{iPar},': ', num2str(PointData(iPar-1))]}]; %#ok<AGROW>
                    end
                otherwise
                    error('Unsupported mode ''%s'' for datatip function.',mode)
            end
        end
    end

% Get numeric value of edit field
    function value = getNum(hObj)
        value = str2double(get(hObj,'String'));
    end

% Set numeric value of edit field
    function setNum(hObj,value,isInteger)
        %         value = num2str(value);
        if nargin<3 || isempty(isInteger)
            isInteger = false;
        end
        
        if isInteger
            set(hObj,'String',sprintf('%i',value));
        else
            set(hObj,'String',sprintf('%.2f',value));
        end
    end

% Callback for edit fields containing floats. Checks if a correct
% number was entered and restricts it to the given bounds.
    function callback_FloatEdit(hObj,~, minVal, maxVal)
        if nargin<3 || isempty(minVal);
            minVal=-inf;
        end
        if nargin<4 || isempty(maxVal);
            maxVal=inf;
        end
        
        % Check if a valid number was entered
        value = str2double(get(hObj, 'String'));
        if isempty(value)
            set(hObj,'ForegroundColor','r');
            set(hObj,'String','INVALID');
            uicontrol(hObj);
        else
            value = max(minVal,value);
            value = min(maxVal,value);
            set(hObj,'ForegroundColor','k');
            set(hObj,'String',sprintf('%.2f',value));
        end
    end

% Callback for edit fields containing integer values. Checks if a correct
% number was entered and restricts it to the given bounds.
    function callback_intEdit(hObj,~, minVal,maxVal)
        if nargin<3 || isempty(minVal);
            minVal=0;
        end
        if nargin<4 || isempty(maxVal);
            maxVal=inf;
        end
        
        value = round(str2double(get(hObj,'String')));
        if isempty(value)
            set(hObj,'ForegroundColor','r');
            set(hObj,'String','INVALID');
            uicontrol(hObj);
        else
            value = max(minVal,value);
            value = min(maxVal,value);
            set(hObj,'ForegroundColor','k');
            set(hObj,'String',sprintf('%i',value));
        end
    end


% The main function of the application. This plays the movie if the timer is running
    function onTimerUpdate(~, ~)
        % Progress frame counter, clip at length of movie and stop at last frame.
        frame = frame+1;
        if(frame >= size(movie,3))
            frame = size(movie,3);
            updateTopText();
            
            stop(h_all.timer);
        end
        set(h_all.slider,'Value',frame);
        updateTopText()
        
        % Skip frame if computer is too slow drawing
        if elapsed_time > timePerFrame
            elapsed_time = elapsed_time - timePerFrame;
            return;
        end
        tic_start = tic;
        updateFrameDisplay();
        elapsed_time = elapsed_time + toc(tic_start)- timePerFrame;
    end

    % Resets the measured elapsed time since last drawn frame when the timer starts
    function onTimerStart(~, ~)
        set(h_all.but_play,'String','Pause');
        elapsed_time = 0;
    end

    function onTimerStop(~, ~)
        set(h_all.but_play,'String','Play');
    end

% Call to display the current frame as selected by the 'frame' variable
% Also this sets and saves the axis states (e.g. for zooming);
    function updateFrameDisplay()
        % Needed to minimize interference with other figures the user
        % brings into focus. It can be that the images are not plotted to
        % the GUI then but to the selected figure window
        set(0,'CurrentFigure',h_main);
        
        % Delete active datatip
        % WARNING this also deletes other hggroup objects associated with the figure which are also invisible and draggable.
        delete(findall(h_main,'Type','hggroup','HandleVisibility','off','Draggable','on'));

        plotFrame(frame);
        
        % Adjust contrast continously if shift key is pressed
        modifiers = get(gcf,'currentModifier');
        shiftIsPressed = ismember('shift',modifiers);
        if(shiftIsPressed)
            autocontrastCallback([],[]);
        end
        
        % Important! Or Matlab will skip drawing entirely for high FPS
        if(MATLAB_2015b_or_newer)
            drawnow nocallbacks; 
        else
            drawnow expose update;
        end
    end

    % Plots contents of the frame with the input index iF
    function plotFrame(iF)
        % Draw frame iF of the movie 
        if imagehandle == -1
            imagehandle = imagesc(movie(:,:,iF)); axis image; colormap gray;
        else
            set(imagehandle,'CData',movie(:,:,iF));
        end
        if use_bw
            colormap gray;
        else
            colormap hot;
        end
        
        % Draw mode dependent data
        switch mode
            case 'movie'
                % Nothing additional to draw.
            case 'candidate'
                if(isempty(candidateData{iF}));
                    if dothandle_cand ~= -1 % Skip uninitialized handles (must be drawn once)
                        set(dothandle_cand,'xdata',[],'ydata',[]);
                    end
                    return % Jump empty frames
                end 
                
                % Plot markers of candidates
                hold on;
                if dothandle_cand == -1 % Draw unitialized handles
                    dothandle_cand = plot(candidateData{iF}(:,1), candidateData{iF}(:,2), 's','Color',marker_color,'MarkerSize',5,'MarkerFaceColor', marker_fill_color');
                else % For initialized handles set their data (MUCH faster than plot)
                    set(dothandle_cand,'xdata',candidateData{iF}(:,1),'ydata',candidateData{iF}(:,2));
                end
                hold off;
                
                
            case 'refinement'
                if(isempty(refinementData{iF}));
                    if dothandle_fit ~= -1 % Skip uninitialized handles (must be drawn once)
                        set(dothandle_fit,'xdata',[],'ydata',[]);
                    end
                    return % Jump empty frames
                end 
                
                % Plot markers of fitted positions
                hold on;
                if dothandle_fit == -1 % Draw unitialized handles
                    dothandle_fit = plot(refinementData{iF}(:,1), refinementData{iF}(:,2), 'o','Color',marker_color,'MarkerSize',5,'MarkerFaceColor', marker_fill_color ,'Linewidth',1);
                else % For initialized handles set their data (MUCH faster than plot)
                    set(dothandle_fit,'xdata',refinementData{iF}(:,1),'ydata',refinementData{iF}(:,2));
                end
                hold off;
                
                
            case 'tracking'
                % Draw the tracks of currently visible particles
                hold on;
                handleNr = 1;
                for iTr = tracksVisibleInFrame{iF}                    
                    % Plot trajectories a) only the last traj_displayLength positions AND  b) up to the current frame
                    mask_toPlot = ((cell_traj{iTr}(:,1)>iF-traj_displayLength) & cell_traj{iTr}(:,1)<=iF);
                    
                    % We use the next free linehandle. If there is no
                    % free handle left, we create a new one on the fly.
                    if (handleNr>numel(linehandles) || linehandles(handleNr) == -1)
                        linehandles(handleNr) = plot(cell_traj{iTr}(mask_toPlot, 2), cell_traj{iTr}(mask_toPlot, 3), '.-','Color',track_colors(iTr,:),'Linewidth',1);
                        linehandleNr_to_TrackNr(handleNr) = iTr;
                        handleNr = handleNr+1;
                    else
                        set(linehandles(handleNr),'xdata',cell_traj{iTr}(mask_toPlot, 2),'ydata', cell_traj{iTr}(mask_toPlot, 3),'Color',track_colors(iTr,:));
                        linehandleNr_to_TrackNr(handleNr) = iTr;
                        handleNr=handleNr+1;
                    end
                end
                % Empty data of unused linehandles
                while(handleNr<= numel(linehandles))
                    if(linehandles(handleNr) ~= -1)                        
                        set(linehandles(handleNr),'xdata',[],'ydata',[]);
                        linehandleNr_to_TrackNr(handleNr) = -1;
                    end
                    handleNr=handleNr+1;
                end
                hold off;
            otherwise
                error('Unkown display mode ''%s''!',mode);
        end
    end

    % Switch play/pause by button
    function playCallback(~, ~)
        if frame == size(movie,3)
            frame = 1;
        end
        if strcmp(get(h_all.timer, 'Running'), 'off')
            start(h_all.timer);
            
            % Disable datatip during playback (otherwise MATLAB tries to
            % put a datatip on the continously updated movie frame, which
            % can lead to errors).
            set(dcm_obj,'Enable','off'); 
        else
            stop(h_all.timer);
            
            % Enable datatip
            set(dcm_obj,'Enable','on'); 
        end
    end

    % Stop playing, adjust contrast manually, continue
    function contrastCallback(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        % We clip the visible display range before the contrast dialog
        % to prevent a warning dialog about display range outside data range
        axes(h_all.axes);
        currImg = movie(:,:,frame);
        zImg = [min(currImg(:)), max(currImg(:))];
        zl = double([max(zImg(1),zl(1)), min(zImg(2),zl(2))]);
        caxis(zl);
        
        % Show contrast dialog, update color axis
        him = imcontrast(h_all.axes);
        uiwait(him);
        zl = caxis;
        
        if isTimerOn
            start(h_all.timer);
        end
    end

    % Stop playing, set contrast to match image min/max values, continue
    function autocontrastCallback(~, ~)
        axes(h_all.axes);
        xl = xlim; % update the axis limits in case the user zoomed
        yl = ylim;
        
        %         currImg = movie(:,:,frame); % Take whole frame for autocontrast
        % Take visible image cutout for autocontrast
        visibleXRange = max(1,floor(xl(1))):min(size(movie,2),ceil(xl(2)));
        visibleYRange = max(1,floor(yl(1))):min(size(movie,1),ceil(yl(2)));
        currImg = movie(visibleYRange,visibleXRange,frame);
                
        selected_method =  get(h_all.popup_autocontrast,'Value');        

        switch selected_method
            case 1 % Focus on spots (upper quartil / 25% of data)
                sortedIntensities = sort(currImg(:));                 
                minIdx = round(0.75*numel(sortedIntensities));
                if minIdx == 0; minIdx = 1; end;
                minval = sortedIntensities( minIdx );
                maxval = sortedIntensities( end );
            case 2 % Adjust contrast to match min/max intensity
                minval = min(currImg(:));
                maxval = max(currImg(:));
            case 3 % 98% range of data (cut upper / lower tails)
                sortedIntensities = sort(currImg(:)); 
                minIdx = round(0.01*numel(sortedIntensities));
                if minIdx == 0; minIdx = 1; end;
                minval = sortedIntensities( minIdx );
                maxval = sortedIntensities( round(0.99*numel(sortedIntensities)));                        
            otherwise
                error('Unknown autocontrast mode!')
        end
        
        zl = double([minval, maxval]);        
        caxis(zl);
    end

    % Plots the distribution/histogram of the parameter selected by popup_distribution
    % Selectable parameters vary depending on data available to the current mode
    function distributionCallback(~, ~)
        selected_parameter = get(h_all.popup_distribution,'Value');
%         choices = get(h_all.popup_distribution,'String');
%         selected_string = choices{selected_parameter};
        
        % distribution to plot, must be synchronized with Popup-Menu (h_all.popup_distribution)
        dataRange = getNum(h_all.edit_distributionRange); % Range of data to histogram
        figure;
        ylabel('frequency');
        switch mode
            case 'candidate'
                % Save concatenated data of all frames (done only once)
                if(isempty(allFramesCandidateData)) 
                    allFramesCandidateData = vertcat(candidateData{:});
                end
                rangedHist(allFramesCandidateData(:,selected_parameter), getNum(h_all.edit_distributionBins),dataRange);
                xlabel(refinementParams{selected_parameter});
            case 'refinement'
                % Save concatenated data of all frames (done only once)
                if(isempty(allFramesRefinementData))
                    allFramesRefinementData = vertcat(refinementData{:});
                end
                rangedHist(allFramesRefinementData(:,selected_parameter), getNum(h_all.edit_distributionBins),dataRange);
                xlabel(refinementParams{selected_parameter});
            case 'tracking'
                rangedHist(trackingData(:,selected_parameter), getNum(h_all.edit_distributionBins),dataRange);
                xlabel(trackingParams{selected_parameter});
            otherwise
                error('Unkown display mode ''%s''!',mode);
        end
    end

    % Switch black-white and hot display mode
    function bwCallback(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        use_bw = ~use_bw;
        
        drawColors(nr_track_colors); % Recompute colors
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay()
        end
        
    end

    % Recompute the colors based on the current background
    function drawColors(num_colors)
        if use_bw
            bg = {'k'}; % background color
        else
            bg = {'r'};
        end
        
        colors = distinguishable_colors(2, bg);
        marker_color = colors(2,:);
        marker_fill_color = colors(1,:);
        if dothandle_cand ~= -1
            set(dothandle_cand,'Color',marker_color, 'MarkerFaceColor',marker_fill_color);
        end
        if dothandle_fit ~= -1
            set(dothandle_fit,'Color',marker_color, 'MarkerFaceColor',marker_fill_color);            
        end
        
        % Draw num_colors colors. If num_colors is less then the number of
        % tracks, we periodically repeat the colors.
        track_colors = repmat( distinguishable_colors(num_colors, bg), ceil(n_tracks/num_colors) ,1);
        track_colors = track_colors(1:n_tracks,:);
    end



    % Update the movie FPS
    function fpsCallback(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        FPS = str2double(get(h_all.edit_FPS, 'String'));
        if isempty(FPS) || FPS<=0
            FPS = 30;
        end
        
        % Timer is limited to 1ms precision
        if FPS>1000
            FPS = 1000;
            warning('Max FPS is 1000 due to timer precision');
        end
        timePerFrame = round(1/FPS*1000)/1000; % limit to millisecond precision
        set(h_all.timer,'Period', timePerFrame);
        set(h_all.edit_FPS,'String',sprintf('%.1f',FPS));
        
        if isTimerOn
            start(h_all.timer);
        end
    end

    % This is called after letting the slider go. We update the frame display
    % once more, otherwise synchronisation issues can occur.
    function sliderCallback(~, ~)
        updateFrameDisplay();
        elapsed_time = 0;
    end

    % This is called continously when dragging the slider
    function updateSlider(~,~)
        % Round slider value, sync text
        frame = round(get(h_all.slider,'Value'));
        set(h_all.slider,'Value',round(frame));
        updateTopText();
        
        % Stop timer
        if strcmp(get(h_all.timer, 'Running'), 'on')
            stop(h_all.timer);
        end
        
        % Skip frame if computer is too slow
        if elapsed_time >timePerFrame
            elapsed_time = elapsed_time - timePerFrame;
            return;
        end
        
        tic_start = tic;
        updateFrameDisplay();
        elapsed_time = elapsed_time + toc(tic_start)- timePerFrame;
    end

    % Sets text on top of the visualizer according to the current frame
    function updateTopText()
        if firstFrame == 1
            set(h_all.toptext,'String',sprintf('Fr. %i/%i',frame,size(movie,3)));
        else
            set(h_all.toptext,'String',sprintf('Fr. %i/%i (%i in movie)',frame,size(movie,3), frame+firstFrame-1));
        end
    end

    % Parse input variables and set them to their default values if they
    % are not given. This is also responsible for loading the movie or a
    % TNT file if these are given as input. Additionally it prepares the
    % trackingData for plotting and disables buttons of modes where the
    % corresponding data is missing.
    function parse_inputs_and_setup(num_argin)
        
        fullPathToThisFile = mfilename('fullpath');
        [TNTpath,~,~] = fileparts(fullPathToThisFile);
        
        % Set unspecified input to empty
        if num_argin<2
            candidateData = [];
        end
        if num_argin<3
            candidateParams = {};
        end
        if num_argin<4
            refinementData = [];
        end
        if num_argin<5
            refinementParams = {};
        end
        if num_argin<6
            trackingData = [];
        end
        if num_argin<7
            trackingParams ={};
        end
        if num_argin<8
            FPS = [];
        end
        if num_argin<9
            is_blocking = [];
        end
        
        if num_argin<10
            firstFrame = [];
        end
        
        % Check if movie, path to movie or path to TNT file was given as first argument
        if(isnumeric(movieOrTNTfile)) % movie was given
            movie = movieOrTNTfile;
        elseif(ischar(movieOrTNTfile)); % path was given
            [~,~,ext] = fileparts(movieOrTNTfile);
            switch ext
                case '.tif' % path to movie was given
                    fprintf('Loading movie from specified path..\n')
                    movie = read_tiff(movieOrTNTfile);
                case '.mat' % path to TNT file was given
                    % Load movie and data
                    fprintf('Loading TNT file..\n')
                    
                    % Before loading, The plugin subfolder path is removed, because loading an XXXoptions
                    % struct with a (main/init/post) function handle where the plugin file exists,
                    % but the subfunction does not (e.g. because it was renamed) throws an
                    % error on loading the file
                    s = warning('off','all');
                        rmpath(genpath([TNTpath,filesep,'plugins']));
                    warning(s);          
                    warning off backtrace
                        TNTdata = load(movieOrTNTfile);
                    warning on backtrace
                    % Add plugin path again (otherwise we loading a TNT
                    % file gives warning that the function to the function
                    % handles could not be found.
                    addpath(genpath([TNTpath,filesep,'plugins']));
                    
                    fprintf('Loading movie specified in TNT file..\n')
                    if(isfield(TNTdata,'firstFrame_lastFrame'))
                        firstFrame = TNTdata.firstFrame_lastFrame(1);
                        movie = read_tiff(TNTdata.filename_movie,false, TNTdata.firstFrame_lastFrame);
                    else
                        movie = read_tiff(TNTdata.filename_movie,false);
                    end
                    transferDatafromTNTdata(TNTdata);
                    
                    % If TNT file is given as first AND second argument, throw error
                    if num_argin>1
                        if(ischar(candidateDataOrTNTfile))
                            error('Invalid input in TNTvisualizer (Two TNT/MAT files).')
                        end
                    end
            end
        end
        
        % Check if candidateData or TNT file was given as second argument
        if num_argin>1
            if ischar(candidateDataOrTNTfile) % TNT file was given
                fprintf('Loading TNT file..\n')
                % Before loading, The plugin subfolder path is removed, because loading an XXXoptions
                % struct with a (main/init/post) function handle where the plugin file exists,
                % but the subfunction does not (e.g. because it was renamed) throws an
                % error on loading the file
                s = warning('off','all');
                    rmpath(genpath([TNTpath,filesep,'plugins']));
                warning(s);
                warning off backtrace
                    TNTdata = load(candidateDataOrTNTfile);
                warning on backtrace
                % Add plugin path again (otherwise we loading a TNT
                % file gives warning that the function to the function
                % handles could not be found.
                addpath(genpath([TNTpath,filesep,'plugins']));
                if(~isequal(size(movie), TNTdata.movieSize))
                   error('Input movie and given TNT file do not fit together! Movie size [%i,%i,%i], TNT file [%i,%i,%i].\n',size(movie),TNTdata.movieSize) 
                end
                if(~isequal(size(movie), TNTdata.movieSize))
                   error('Input movie and given TNT file do not fit together! Movie size [%i,%i,%i], TNT file [%i,%i,%i].\n',size(movie),TNTdata.movieSize) 
                end                
                transferDatafromTNTdata(TNTdata);
            elseif iscell(candidateDataOrTNTfile)
                candidateData = candidateDataOrTNTfile;
            end            
        end
        
        
        % -- Initialize variables left empty --
        % Is candidate data available?
        if isempty(candidateData)
            set(h_all.button_candidateMode,'Enable','off');
        else
            candidateParams = checkParameterDescription(candidateData, candidateParams);
            mode = 'candidate';
        end
        
        if isempty(candidateParams)
            candidateParams = {};
        end
        
        % Is refinement data available?
        if isempty(refinementData)
            set(h_all.button_refinementMode,'Enable','off');
        else
            refinementParams = checkParameterDescription(refinementData, refinementParams);
            mode = 'refinement';
        end
        
        if isempty(refinementParams)
            refinementParams = {};
        end
        
        % Is tracking data available?
        if isempty(trackingData)
            set(h_all.button_trackingMode,'Enable','off');
            
            id_tracks = []; % Note: (in case track IDs go from 1 to N without missing numbers, the track index is identical to the tracks ID.
            n_tracks = 0;
        else
            trackingParams = checkParameterDescription(trackingData, trackingParams);
            mode = 'tracking';
            
            id_tracks = unique(trackingData(:,1)); % Note: (in case track IDs go from 1 to N without missing numbers, the track index is identical to the tracks ID.
            n_tracks = numel(id_tracks);
        end
        % Prepare tracking data for use with the visualizer
        % Put data for each track into its own cell. This is much more efficient comfortable for plotting.
        cell_traj = cell(n_tracks ,1);
        cnt = 1;
        for iTrack = 1:n_tracks
            cell_traj{iTrack} = trackingData( trackingData(:,1)== id_tracks(cnt) , 2:end);
            cnt = cnt+1;
        end
        
        
        if isempty(trackingParams)
            trackingParams ={};
        end
        
        if isempty(FPS)
            FPS = 30;
        end
        
        if isempty(is_blocking)
            is_blocking = false;
        end
        
        if isempty(firstFrame)
            firstFrame = 1;
        end
    end

    % Transfer data from TNTfile into our local GUI variables
    function transferDatafromTNTdata(TNTdata)
        if(isfield(TNTdata,'candidateData'))
            candidateData = TNTdata.candidateData;
            candidateParams = TNTdata.candidateOptions.outParamDescription;
        else
            candidateData = {};
            candidateParams = {};
        end
        if(isfield(TNTdata,'refinementData'))
            refinementData = TNTdata.refinementData;
            refinementParams = TNTdata.refinementOptions.outParamDescription;
        else
            refinementData = {};
        end
        if(isfield(TNTdata,'trackingData'))
            trackingData = TNTdata.trackingData;
            trackingParams = TNTdata.trackingOptions.outParamDescription;
        else
            trackingData = {};
            trackingParams = {};
        end
        if(isfield(TNTdata,'firstFrame_lastFrame'))
            firstFrame = TNTdata.firstFrame_lastFrame(1);
        end
    end

    % Check if the number of parameters in the data and their description match.
    % If there are not enough descriptions, pad with '<Unknown>', if there are
    % too much, truncate.
    function description = checkParameterDescription(data, description)
        % Find number of parameters in fitData
        if iscell(data)
            nrParams = 0;
            for iFrame = 1:size(data,1)
                if ~isempty(data{iFrame})
                    nrParams = size(data{iFrame},2);
                    break;
                end
            end
        else
            nrParams = size(data,2);
        end
        
        % If there is no data, do nothing
        if nrParams == 0
            return
        end
        
        % Make parameter description the same size as the number of params
        if isempty(description)
            description = repmat({'<Unknown>'}, nrParams,1);
        else
            if numel(description) > nrParams
                description = description(1:nrParams);
            elseif numel(description) < nrParams
                tmp = description;
                description = cell(nrParams,1);
                description(:) = {'<Unknown>'};
                description(1:numel(tmp)) = tmp(1:numel(tmp));
            end
        end
        
    end

    % Cleanup function. This is neccessary to delete the timer!
    function onAppClose(~, ~)
        if strcmp(get(h_all.timer, 'Running'), 'on')
            stop(h_all.timer);
        end
        delete(h_all.timer);
        delete(h_main);
    end


%%  Tracking mode related functions

    % Compute which tracks are visible in each frame
    % This is done to iterate only over visible tracks when plotting. Otherwise
    % a large number of tracks in the movie slows us down, even if we skip them
    % on the fly (because we have to check each track for every frame)
    %
    % Sets:
    %  tracksVisibleInFrame
    %  maxNr_visibleTracksInFrame
    function compute_tracksVisibleInFrame()
        nr_frames = size(movie,3);
        tracksVisibleInFrame = cell(nr_frames,1);
        tracksVisibleInFrame_bool = false(n_tracks,nr_frames);
        
        %Note: We need the bool to compute this fast (at the cost of memory)
        % but convert the data to a cell array for using it later
        for iTrack = 1:n_tracks
            track_firstFrame = cell_traj{iTrack}(1,1);
            track_lastFrame = cell_traj{iTrack}(end,1);

            tracksVisibleInFrame_bool(iTrack, track_firstFrame:min(track_lastFrame+traj_lifetime, nr_frames)) = true;
        end
        
        % Convert bool matrix to indices
        for iFrame = 1:nr_frames
            tracksVisibleInFrame{iFrame} = find(tracksVisibleInFrame_bool(:,iFrame)).';
        end

        [~,nr_VisibleTracksInFrame] = cellfun(@size,tracksVisibleInFrame);
        maxNr_visibleTracksInFrame = max(nr_VisibleTracksInFrame);
    end

    % Simply plots the frame with the maximum number of visible tracks,
    % which initializes all linehandles needed for displaying the movie. If
    % this is not called, linehandles are simply initialized on the fly
    % (which might cause jerky playback on playing the movie)
    %
    % Make sure updateFrameDisplay() is called afterwards, so that the
    % current frame is displayed again!
    function initializeLinehandles()
       [~, maxFrame] = max(nr_VisibleTracksInFrame);
       plotFrame(maxFrame);
    end

     % Sets the traj_lifetime variable if it changes in the GUI
    function callback_TrajLifetime(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        traj_lifetime = round(str2double(get(h_all.edit_lifetime,'String')));
        if traj_lifetime<=0 || isempty(traj_lifetime)
            traj_lifetime = 0;
        end
        set(h_all.edit_lifetime,'String',sprintf('%i%',traj_lifetime));
        
        % We have to update which tracks are visible in which frames
        % This alters the number of linehandles which are needed for
        % display (the maximum number of concurrently visible tracks)
        compute_tracksVisibleInFrame();
        resetGraphics();
        initializeLinehandles();
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end

    % Sets the traj_displayLength variable if it changes in the GUI
    function callback_dispLength(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        traj_displayLength = round(str2double(get(h_all.edit_trajDisplayLength,'String')));
        if traj_displayLength<=0 || isempty(traj_displayLength)
            traj_displayLength = 0;
        end
        set(h_all.edit_trajDisplayLength,'String',sprintf('%i%',traj_displayLength));
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end

    % The user entered a different number of colors -> update color pool
    function callback_trackColors(~, ~)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        nr_track_colors = round(str2double(get(h_all.edit_colors,'String')));
        if nr_track_colors<=0 || isempty(nr_track_colors)
            nr_track_colors = 1;
        end
        set(h_all.edit_colors,'String',sprintf('%i%',nr_track_colors));
        
        drawColors(nr_track_colors);
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end
end


%% --- General functions ---

% Adds pathes needed for the visualizer.
function addPathsVisualizer()
    fullPathToThisFile = mfilename('fullpath');
    [path,~,~] = fileparts(fullPathToThisFile);
    addpath(genpath([path,filesep,'subfun']));
    addpath(genpath([path,filesep,'external']));
    addpath(genpath([path,filesep,'helper']));
end

% Function that cuts data from upper and lower tails of the distribution
% keeping at least 'percentOfData' percent of all values.
function [freq, centers] = rangedHist(data, nbins, percentOfData)
if(percentOfData<100)
    % Limits for the cumulative density function (which goes from 0 to 1)
    lower_limit = (1-percentOfData/100)/2; % below this we throw away
    upper_limit = 1-(1-percentOfData/100)/2; % above this we throw away
    
    data = sort(data);
    minIdx = round(lower_limit*numel(data));
    if minIdx < 1
       minIdx = 1; 
    end
    maxIdx = round(upper_limit*numel(data));
    if maxIdx < 1
       maxIdx = 1; 
    end
    
    data = data(minIdx:maxIdx);
end

freq = [];
centers = [];
switch(nargout)
    case 0
        hist(data,nbins);
    case 1
        freq = hist(data,nbins);
    case 2
        [freq,centers] = hist(data,nbins);
end

end
