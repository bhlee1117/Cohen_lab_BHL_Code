% example analysis using vm class
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%

% open latest movie recorded in r:\
foldername = 'R:';

mdir = dir(fullfile(foldername,'**'));
mdir = mdir(cell2mat({mdir.isdir}));
mdir = mdir(arrayfun(@(x)~isequal(x.name,'.'),mdir));
mdir = mdir(arrayfun(@(x)~isequal(x.name,'..'),mdir));
mdir = mdir(arrayfun(@(x)~isequal(x.name,'temp'),mdir));
mdir = mdir(arrayfun(@(x)~isequal(x.name,'rec'),mdir));
mdir = mdir(arrayfun(@(x)~isequal(x.name,'Temp'),mdir));
[d,ix] = sort(arrayfun(@(x)x.name,mdir,'uni',false));
mdir = mdir(ix(end));

binname = fullfile(foldername,mdir.name,'Sq_camera.bin');
dcimgname = fullfile(foldername,mdir.name,'Sq_camera.dcimg');
for it = 1:600
    if ~exist(binname,'file')
        pause(.1)
    else
        break
    end
end
for it = 1:600
    if exist(dcimgname,'file')
        pause(.1)
    else
        break
    end
end

disp(mdir.name)
clear mov
mov = vm(fullfile(foldername,mdir.name))

useframes = 3:62;
hadtraces = hadamard_bincode_nopermutation(59)'*2-1;
moviefixsc(mov(useframes)*hadtraces)

% [rois, intens] = nested_clicky_rects(mov');
%%

figure windowstyle docked
intens = apply_clicky_faster(mov',rois);
title(strrep(mdir.name,'_','\_'))
saveas(gcf,fullfile(foldername,mdir.name,[datestr(now,'HHMMSS') '_clicky.fig']))
drawnow

figure windowstyle docked
semilogy(abs(fft(intens)))
title(strrep(mdir.name,'_','\_'))
saveas(gcf,fullfile(foldername,mdir.name,[datestr(now,'HHMMSS') '_auto_spectrum.fig']))
ylim(10.^[-1 7])
drawnow
%%

% useframes = 3:62;
% hadtraces = hadamard_bincode_nopermutation(59)'*2-1;
% moviefixsc(mov(useframes)*hadtraces)
movmean = mov.mean;
saveastiff(uint16(movmean),fullfile(foldername,mdir.name,'mean.tif'));
figure
moviefixsc(mov-movmean)
title(mdir.name,'Interpreter','none')
figure
plot(mov.frameAverage)
title(mdir.name,'Interpreter','none')
saveas(gcf,fullfile(foldername,mdir.name,'frameAverage.fig'));
saveas(gcf,fullfile(foldername,mdir.name,'frameAverage.png'));

%%
figure windowstyle docked
plot(mov.frameAverage)
% moviefixsc(mov([1:144 (-143:-0)+end]).evnfun(@mean,12)*hadtraces)
title(mdir.name)
useframes = [2+(1:132) 166+(1:132)]'+300*(0:5);

% clf
% plot(mov.frameAverage)
% hold on
% for it = 1:size(useframes,2)
%     plot(useframes(:,it),mov(useframes(:,it)).frameAverage)
% end
%%
% plot(mov(useframes).frameAverage)
moviefixsc(mov(useframes)*hadtraces)
