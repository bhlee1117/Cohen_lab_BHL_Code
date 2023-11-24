function testInPlaceSimple

% Copyright 2007 The MathWorks, Inc.

% INTERNAL USE ONLY
%  Using "feature mtic, mtoc"
%  Alternatively, bring up Task Manager

%% Create Variables (Pre-Allocate)
disp('Creating Variables');drawnow;
x = rand(5000);
y = rand(5000);
pause(3);

%% In-Place
disp('in-place');drawnow;
feature mtic;
x = 2*x + 3;
feature mtoc
pause(3);

%% NOT In-Place
disp('not in-place');drawnow;
feature mtic
y = 2*x + 3;
feature mtoc
pause(3);