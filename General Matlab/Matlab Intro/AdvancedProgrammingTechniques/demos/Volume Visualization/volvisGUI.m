function volvisGUI(x,y,z,v)

% Copyright 2007 The MathWorks, Inc.

% volvisGUI provides interactive volume visualization
%
% Ex:
% [x,y,z,v] = flow;
% volvisGUI(x,y,z,v)

%% Default Input
if nargin == 0
  [x,y,z,v] = flow;
end

%% Initalize visualization
figure;
fh = volumeVisualization(x,y,z,v);
fh.addSlicePlane(fh.xMin);

%% Add uicontrol
hSlider = uicontrol(...
  'Units'   , 'normalized'          , ...
  'Position', [.75 .05 .2 .05]      , ...
  'Style'   , 'slider'              , ...
  'Min'     , fh.xMin               , ...
  'Max'     , fh.xMax               , ...
  'Value'   , fh.xMin               , ...
  'Callback', @updateSliderPosition );

  % Slider Callback Function
  function updateSliderPosition(varargin)
    fh.deleteLastSlicePlane();
    x = get(hSlider,'Value');
    fh.addSlicePlane(x);
  end

end
