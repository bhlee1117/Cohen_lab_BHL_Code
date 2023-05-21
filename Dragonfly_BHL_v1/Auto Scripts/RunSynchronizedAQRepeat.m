
nrpt = 20;
tpause = 1;
%%
for ii=1:nrpt
    hDragonflyApp.RunSynchronizedAQ([])
    pause(tpause)
end