function voltage_sub=get_subthreshold(voltage,spike,dialate_size)
if nargin<3
    dialate_size=3;
end

se = strel('line',dialate_size,0);
spike_dialated=imdilate(spike>0,se);
voltage(spike_dialated==1)=NaN;

voltage_sub=zeros(size(voltage));
for i=1:size(voltage,1)
    t=find(~isnan(voltage(i,:)));
    voltage_sub(i,:)=interp1(t,voltage(i,t),[1:size(voltage,2)],'linear');
end
end