function trace_blueoff=get_blueoffTrace(trace,Blue,erode_sz)
if length(trace)==length(Blue)

t=[1:length(trace)];
blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,erode_sz), zeros(1, erode_sz)]);
trace(~blueOff2)=NaN;
trace_lp=movmean(trace,150,1,'omitnan');
trace_blueoff=interp1(t(blueOff2),trace_lp(blueOff2),t,'linear')';

else
    error('Length of trace and blue input is different')
end