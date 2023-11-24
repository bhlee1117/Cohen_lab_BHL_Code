%% Simple SPMD Example
% This example illustrates our spmd capabilites. 

%% open up a matlabpool

if matlabpool('size') == 0
    matlabpool open 2
else
    warning('MATLAB:poolOpen', 'matlabpool already open')
end

%% serial part of code

n = 100;

%% simple spmd block

spmd
    a = rand(n,n);
    display(size(a))
    display(a(1:2,1:2))
end

%% spmd with distributed arrays (creating)

spmd
    a = rand(n,n,codistributor); % uses default distribution
    display(size(getLocalPart(a))); 
end

%% spmd with distributed arrays (duplicated)

spmd
    b = ones(n,n); 
    a = codistributed(b) ;% uses default distribution
    display(size(b));
    display(size(a)); 
end

%% spmd with distributed arrays (variant)

spmd
    d = rand(n,n); 
    codist = codistributor1d(0,[],[n, numlabs*n]) ;% uses default distribution
    e = codistributed.build(d,codist);
    display(size(d));
    display(size(e)); 
end

%% Data stays on workers between spmd blocks

spmd
    d = svd(a);
    display(max(d))
end

%% Can access the data by referencing into it, now copies to client

display(d{1})

%% Close the matlabpool 

if matlabpool('size') ~= 0
    matlabpool close
else
    warning('MATLAB:poolClose', 'matlabpool already closed')
end


