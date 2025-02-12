function mov=readBinMov_BHL_multiple(filepath,n,time_read,time_segment)

if nargin<2
    n=3;
end
load([filepath '/output_data.mat'])
sz=double(Device_Data{n}.ROI([2 4]));

readall=0;
if nargin<3
    time_read=[1:max(Device_Data{1, 2}.Counter_Inputs.data)];
    readall=1;
    time_segment=15000;
end

isMCtrace=dir(fullfile(filepath,['mc_*.bin']));
mov=[];
file_read=ceil(time_read/time_segment);
frame_read=mod(time_read,time_segment);
frame_read(frame_read==0)=time_segment;

if ~isempty(isMCtrace)

    for j=unique(file_read)
    disp(['Loading ' fullfile(filepath,isMCtrace(j).name)])
mov=cat(3,mov,double(readBinMov_times(fullfile(filepath,isMCtrace(1).name),sz(2),sz(1),frame_read(file_read==j))));    
    end
else
mov=double(readBinMov_times(fullfile(filepath,'frames1.bin'),sz(2),sz(1),time_read));
end

end