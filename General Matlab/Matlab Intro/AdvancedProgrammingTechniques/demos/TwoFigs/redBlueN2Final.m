function redBlueN2Final

% Copyright 2007 The MathWorks, Inc.

h(1) = figure('Color', 'k');movegui northwest
h(2) = figure('Color', 'k');movegui northeast
red  = uicontrol(               ...
  'String'  , 'red'           , ...
  'Position', [20 40 60 20]   , ...
  'Parent'  , h(1)            , ...
  'Callback', @redButtonPress );
blue = uicontrol(               ...
  'String'  , 'blue'          , ...
  'Position', [20 20 60 20]   , ...
  'Parent'  , h(2)            , ...
  'Callback', @blueButtonPress);

  function redButtonPress(varargin)
    % varargin input required in callbacks even though the info isn't used,
    % because matlab will issue the callback with at least 3 inputs:
    %         the handle to the object triggering the callback,
    %         a dummy variable reserved for future use, and
    %         event information
    % with this nested function example, we don't need to use that information
    % but we need to be able to accept the inputs
    set(red , 'Enable', 'off');
    set(blue, 'Enable', 'on' );
    set(h   , 'Color' , 'r'  );

  end

  function blueButtonPress(varargin)

    set(red , 'Enable', 'on' );
    set(blue, 'Enable', 'off');
    set(h   , 'Color' , 'b'  );

  end

end
