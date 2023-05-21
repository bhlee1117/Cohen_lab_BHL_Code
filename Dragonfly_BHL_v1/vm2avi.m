function vm2avi(obj,fname,framerate)
% function vm2avi(obj,fname,framerate)
% save vm to avi, no compression
if nargin > 1
    pname = '';
else
    [fname, pname] = uiputfile('*.avi');
end
fullname = fullfile(pname,fname);
tic
fprintf('saving %s ...\n',fullname)
v = VideoWriter(fullname);
v.FrameRate = framerate;
v.Quality = 95;
v.open

dat = single(permute(obj.toimg.data,[1 2 4 3]));
dat = dat/max(dat(:));
writeVideo(v,dat)
v.close
disp(['saving took ' num2str(toc) ' s']);

end