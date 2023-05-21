ROI=[];
intensities=[];

[fpath] = uigetdir; if fpath == 0, return;end
cd(fpath); list=dir;
for i=1:length(list)
    if  list(i).isdir && ...
        (~(strcmp(list(i).name,'.')||strcmp(list(i).name,'..')))
    
        if isempty(ROI)
            mov = vm(fullfile(list(i).folder,list(i).name));
            mov = mov(1:80);
            % mov=setSaturationLimits(mov,[0 700]);
            figure;moviesc(mov)%imshow(mov.mean,[])

            [ROI,intens]=clicky_autobg(mov.data,mov(1:end/2).mean);
            intensities = [intensities intens];
        else
            mov = vm(fullfile(list(i).folder,list(i).name));

            intens = apply_clicky_autobg(mov,ROI,'no');
            intensities = [intensities intens];
        end
    end 
end

offsets = range(intensities,1);offsets(end) = 0;offsets = circshift(offsets,1);
% figure;plot(intensities+offsets*triu(ones(size(intensities,2))))
t = 0:.03:(length(intensities)-1)*.03;
figure;plot(t,intensities)
xlabel('Time (s)')
% set(gca,'ColorOrder',jet(size(intensities,2)))
% ylabel('{\Delta}F/F')
figure;plot(intensities)