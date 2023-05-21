function dragonfly_shutdown(app,event)

%             cameraMex shutdown
%-------- move snap logs ---------------
if exist('Snap_logs','dir')

    try 
        fname = dir('log*'); fname = fname.name;
        movefile (fname,['Snap_logs/' fname]);
    end
else
    mkdir('Snap_logs')
    try
    fname = dir('log*'); fname = fname.name;
     movefile (fname,['Snap_logs/' fname]);
    end
end
%-------- move snap logs END---------------

try %#okay no catch
    app.camApp.ShutdownButton.ButtonPushedFcn(app,event); 
end

try %#okay no catch
    app.dmdApp.UIFigureCloseRequest(app.dmdApp,event); 
end

try %#okay no catch
    app.slmApp.UIFigureCloseRequest(app.slmApp,event); 
end


try %#okay no catch
    app.LimeLaserstatusButton.Value = false;
    app.LimeLaserstatusButtonValueChanged;
    
    app.Blue488nmButton.Value = false;
    app.Blue488nmButtonValueChanged;
    
    app.ODWheel.moveTo('0');
end

try
    app.motorW.switchTo('RCaMP')
end
evalin( 'base', 'clear hCameraApp' )
evalin( 'base', 'hDaq.syncSession.release' )
evalin( 'base', 'clear hDaq' )
evalin( 'base', 'clear hLimeLaser' )
evalin( 'base', 'clear hAOTF' )
evalin( 'base', 'clear hODWheel' ) 
try %#okay no catch
    evalin( 'base', 'close(hZMotorFig)' );
    evalin( 'base', 'close(hWMotorFig)' );
end
evalin( 'base', 'clear hZMotorFig' )
evalin( 'base', 'clear hWMotorFig' )
evalin( 'base', 'clear hMotorZ' )
evalin( 'base', 'clear hMotorW' )
try %#okay no catch
    app.closeRequestFcn(event); 
end 
evalin( 'base', 'clear hDragonflyApp' )          