function redBlueN

% Copyright 2007 The MathWorks, Inc.

h    = figure('Color', 'k');
red  = uicontrol(               ...
  'String'  , 'red'           , ...
  'Position', [20 40 60 20]   , ...
  'Callback', @redButtonPress );
blue = uicontrol(               ...
  'String'  , 'blue'          , ...
  'Position', [20 20 60 20]   , ...
  'Callback', @blueButtonPress);


  function redButtonPress(varargin)

    set(red, 'Enable', 'off');
    set(blue, 'Enable', 'on');
    set(h, 'Color', 'r');
  end

  function blueButtonPress(varargin)

    set(red, 'Enable', 'on');
    set(blue, 'Enable', 'off');
    set(h, 'Color', 'b');
  end

end
