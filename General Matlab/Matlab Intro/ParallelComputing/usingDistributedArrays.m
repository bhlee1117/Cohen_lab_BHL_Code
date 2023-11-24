%% An introduction to creating and using distributed arrays
%% Opening a MATLAB pool
if ~matlabpool('size')
    matlabpool open 4
end

%% Creating an array

a = rand(10000,10000);

%% Creating a distributed array

D = distributed.rand(10000);  % Data is created and stored on the workers.
D2 = distributed.rand(1000);  % Data is created and stored on the workers.

%% Checking the sizes and classes
% Both D and D2 are the same size in the client workspace, since data is
% held on the workers.
whos

%% Doing math on distributed arrays

b = distributed.rand(10000, 10000); % Created on the workers
x = D + b;

%% Seeing what you can do with distributed arrays

methods(distributed)  

%% Plotting distributed arrays
% The "plot" command is overloaded to work with distributed arrays.
% This brings data back to the client.

plot(x(1:1000:end,1:1000:end))  % brings back every 100th element

%% Gathering the data to a client-side double array
% Note that Xgather is much larger than x in the client workspace.
Xgather = gather(x(1:1000:end,1:1000:end));
whos('Xgather', 'x')

%% Closing MATLAB pool
matlabpool close