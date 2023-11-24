%% Memory Overhead for Different Kinds of Arrays
% How much memory do different datatypes consume?

% Copyright 2007 The MathWorks, Inc.

%% Setup 
clear all
clc

%%
d = [1 2]
dcell = {d}
dstruct.d = d

%% Look at whos Output
whos
