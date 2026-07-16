function nfr = ms2frames(ms, dt_s)
% MS2FRAMES  Convert a duration in milliseconds to a number of frames.
%   nfr = ms2frames(ms, dt_s)
%   dt_s : seconds per frame (e.g. median(diff(t_ax)), with t_ax in seconds).
nfr = max(0, round(ms/1000/dt_s));
end
