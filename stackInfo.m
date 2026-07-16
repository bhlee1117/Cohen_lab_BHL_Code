function I = stackInfo(S, offsetVal)
% STACKINFO  Pack a stacker's reload metadata (see loadSTAevent).
I.roi=S.roi; I.cls=clsName(S.cls);
I.win=S.win; I.winLen=S.win(1)+S.win(2)+1;
I.H=S.H; I.W=S.W; I.pixels=S.H*S.W;
I.offset=offsetVal;
I.centerFrames=S.centerFrames;   % global movie frame of each stored window (stored order)
I.binFiles=S.binFiles;
I.evPerFile=S.evPerFile;
I.n=S.cnt;
end
