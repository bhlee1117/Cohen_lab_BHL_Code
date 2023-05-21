function ringBell()
% 
% Plays a short audio clip of a bell ringing. 
[Y,Fs]=audioread(fullfile('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Code/','ringBell.ogg'));
sound(Y,100000);
end