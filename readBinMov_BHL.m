function mov=readBinMov_BHL(filepath,n,time_read)

if nargin<2
    n=3;
end
load([filepath '/output_data.mat'])
sz=double(Device_Data{n}.ROI([2 4]));

if nargin<3
    time_read=[1:max(Device_Data{1, 2}.Counter_Inputs.data)];
end


isMCtrace=dir(fullfile(filepath,['mc_*.bin']));
if ~isempty(isMCtrace)
    disp(['Loading ' fullfile(filepath,isMCtrace(1).name)])
mov=double(readBinMov_times(fullfile(filepath,isMCtrace(1).name),sz(2),sz(1),time_read));    
else
mov=double(readBinMov_times(fullfile(filepath,'frames1.bin'),sz(2),sz(1),time_read));
end

end