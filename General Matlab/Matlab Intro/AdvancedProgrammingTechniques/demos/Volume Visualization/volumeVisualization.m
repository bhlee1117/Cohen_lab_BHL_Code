function fh = volumeVisualization(x,y,z,v)

% Copyright 2007 The MathWorks, Inc.

%% Data Initialization
% Create a variable to store pointers to the various planes
hAxis = [];         %initialize handle to axis
hSlicePlanes = [];  %initialize handle to slice plane

% Create data for generic slice through yz-plane
[yd,zd] = meshgrid(linspace(min(y(:)),max(y(:)),100), ...
  linspace(min(z(:)),max(z(:)),100));

%% Plot the volume initially
initDisplay()

%% Nested Functions

  % Add a slice plane at specified x-coordinate.
  function addSlicePlane(xLoc)
    xd            = xLoc*ones(size(yd));
    newSlicePlane = slice(hAxis, x, y, z, v, xd, yd, zd);
    hSlicePlanes  = [ hSlicePlanes, newSlicePlane ];
    set(newSlicePlane            , ...
      'FaceColor'      , 'interp', ...
      'EdgeColor'      , 'none'  , ...
      'DiffuseStrength', .8       );
  end

  % Delete the last slice plane
  function deleteLastSlicePlane()
    if ~isempty(hSlicePlanes)
      delete(hSlicePlanes(end));
      hSlicePlanes = hSlicePlanes(1:end-1);
    end
  end

  % Initialize Display
  function initDisplay()
    % initDisplay  Initialize Display

    % Draw back and bottom walls
    if isempty(hAxis) || ~ishandle(hAxis)
      hAxis = gca;
      hold on;
    end
    hx = slice(hAxis, x, y, z, v,max(x(:)),       [],       []) ;
    hy = slice(hAxis, x, y, z, v,       [],max(y(:)),       []) ;
    hz = slice(hAxis, x, y, z, v,       [],       [],min(z(:))) ;

    % Make everything look nice
    set([hx hy hz], 'FaceColor', 'interp', 'EdgeColor' , 'none')
    set(hAxis     , 'FontSize' , 18      , 'FontWeight', 'Bold');
    xlabel('X');ylabel('Y');zlabel('Z')
    daspect([1,1,1])
    axis tight
    box on
    view(-38.5, 16)
    colormap(jet(128))
  end

%% Export structure
fh.addSlicePlane        = @addSlicePlane       ;
fh.deleteLastSlicePlane = @deleteLastSlicePlane;
fh.xMin                 = min(x(:))            ;
fh.xMax                 = max(x(:))            ;

end