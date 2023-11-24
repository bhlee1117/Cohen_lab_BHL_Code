%% Simple MPI Example
% Sending Data from Lab1 to Lab2 using LabSend
% Lab2 Receives Data from Lab1 using LabReceive
%
% This is running in an spmd block but the code inside the blcok
% could run in a paralleljob or in pmode. 

clear all

source = 1;
destination = 2;

%% Open up a matlabpool

if matlabpool('size') == 0
    matlabpool open 2
else
    warning('MATLAB:poolOpen', 'matlabpool already open')
end

%%  Running mpi commands in spmd block
% The code inside the spmd block would also run in a paralleljob or in pmode. 

N = 10;

spmd
    if labindex == source
        % Send some data from source lab
        testData = rand(N,N);
        labSend(testData, destination); 
        
    elseif labindex == destination
        % Receive on destination lab
        recvdata = labReceive(source);
        disp(recvdata)
        
    end
end

%% Access data on the client 

disp(recvdata)


%% Close the matlabpool 

if matlabpool('size') ~= 0
    matlabpool close
else
    warning('MATLAB:poolClose', 'matlabpool already closed')
end

