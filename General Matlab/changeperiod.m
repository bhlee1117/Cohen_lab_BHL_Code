function out = changeperiod(t1, y1, p1, p2);
% function out = changeperiod(t1, y1, p1, p2);
% Original data is y1, sampled at times t1.
% output data is ???
% returns data sampled at times t1, with period scaled by p2/p1.

t2 = t1*p2/p1;
out = interp1(t2, y1, t1, 'linear', 0);
