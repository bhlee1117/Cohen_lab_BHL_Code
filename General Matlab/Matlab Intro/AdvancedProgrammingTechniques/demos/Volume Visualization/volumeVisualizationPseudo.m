function fh = volumeVisualizationPseudo(x,y,z,v)

% Copyright 2007 The MathWorks, Inc.

% volumeVisualizationPseudo   Pseudo-code for volumeVisualization

%% Initialize
initDisplay(x,y,z,v)

%% Nested Functions

  % Add slice plane
  function addSlicePlane(xLoc)
  end

  % Delete Slice plane
  function deleteLastSlicePlane()
  end

  % Initialize Display
  function initDisplay(x,y,z,v)
  end

% Export structure
fh.addSlicePlane        = @addSlicePlane    ;
fh.deleteLastSlicePlane = @deleteLastSlicePlane;
fh.xMin                 = min(x(:))        ;
fh.xMax                 = max(x(:))        ;
end
