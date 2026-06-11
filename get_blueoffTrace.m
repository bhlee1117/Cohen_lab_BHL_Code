function [trace_blueoff blueOff2 blueOff3]=get_blueoffTrace(trace,Blue,erode_sz,erode_sz2)
% GET_BLUEOFFTRACE  Extract and interpolate a trace during blue-light OFF periods.
%
% Usage:
%   [trace_blueoff, blueOff2, blueOff3] = get_blueoffTrace(trace, Blue, erode_sz)
%   [trace_blueoff, blueOff2, blueOff3] = get_blueoffTrace(trace, Blue, erode_sz, erode_sz2)
%
% Description:
%   Identifies time points where the blue light is OFF (Blue == 0), then
%   erodes the edges of each OFF epoch to exclude transient artifacts near
%   blue-light transitions. The valid (eroded) blue-OFF samples are smoothed
%   with a moving average and then linearly interpolated across the full
%   time vector to produce a continuous baseline trace.
%
% Inputs:
%   trace      - Nx1 signal trace
%   Blue       - Nx1 blue-light signal (0 = off, >0 = on)
%   erode_sz   - Number of samples to erode from the front of each blue-OFF epoch
%                (removes artifact immediately after blue turns off)
%   erode_sz2  - Number of samples to erode from the back of each blue-OFF epoch
%                (removes artifact immediately before blue turns on; default = 0)
%
% Outputs:
%   trace_blueoff - Nx1 interpolated trace using smoothed blue-OFF samples
%   blueOff2      - Nx1 logical mask: blue-OFF epochs with front edge eroded
%   blueOff3      - Nx1 logical mask: blue-OFF epochs with both edges eroded
if nargin<4
    erode_sz2=0;
end
if length(trace)==length(Blue)

t = 1:length(trace);
blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,erode_sz), zeros(1, erode_sz)]);
blueOff3 = imerode(blueOff2, [zeros(1, erode_sz2), ones(1,erode_sz2)]);
trace(~blueOff3) = NaN;
trace_lp = movmean(trace, 150, 1, 'omitnan');
trace_blueoff = interp1(t(blueOff3), trace_lp(blueOff3), t, 'linear')';

else
    error('Length of trace and blue input is different')
end