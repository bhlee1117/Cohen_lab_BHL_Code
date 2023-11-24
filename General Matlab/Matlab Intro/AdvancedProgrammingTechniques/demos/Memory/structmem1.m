%% Example: Managing Structure Memory 

% Copyright 2007 The MathWorks, Inc.

clear all, close all

%% Create a structure
% (4000*3200 elts) * (8 bytes/elt) = 100,000 KB ~= 100MB
% Could be interesting to bring up task manager here
s.A = rand(4000,3200);

%%
s.B = rand(4000,3200);

%% Copy the structure
% What will happen to memory?
% Watch task manager here
sNew = s;       

%% Modify the structure
% Can you explain what happened to memory?
s.A(1,1) = 17;   

%%
sNew.B(1,1) = 0;