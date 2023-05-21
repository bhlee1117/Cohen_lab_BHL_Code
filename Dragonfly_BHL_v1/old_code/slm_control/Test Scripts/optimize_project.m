xyz_offset = [0 0 0];
% xyz_offset = best_xyz_offset;
% xyz_offset(:,3) = 5
xyz = slmlocs - xyz_offset;
ttd = [
    0 1 1/3 % 0 sind(60) cosd(60)
    1 0 0
    0 -1/20 1 % cosd(60)/100 sind(60)
    ]*xyz';
approxSLM = gs(commonz, zpolys, ttd, 100, false, true, false, true, true);%, (amps./max(amps)).^2);

slmph1 = angle(approxSLM);
%         slmph1 = (sign(slmph1)/2+.5)*128;
slmph1 = slmph1*128/pi+128;
slmph1 = max(min(slmph1,255),0);
slmph1 = uint8(slmph1);
figure windowstyle docked
imshow(slmph1,[])
colorbar
projected = true;
prefix = [datestr(now,'HHMMSS') '_optimize_project'];
% a = dbstack;
% if ~isempty(dbstack)
if exist('best_xyz_offset','var')
    save(fullfile(ramfovpath,[prefix '.mat']),'slmlocs','slmph1','xyz','ttd','lastClickedX','lastClickedY','priorClickedX','priorClickedY','best_xyz_offset')
else
    save(fullfile(ramfovpath,[prefix '.mat']),'slmlocs','slmph1','xyz','ttd','lastClickedX','lastClickedY','priorClickedX','priorClickedY')
end
% calllib('BlinkHdmiSdk', 'Write_image', sdk, slmph1, 1);
fullscreen(rot90(slmph1),3)
%%

