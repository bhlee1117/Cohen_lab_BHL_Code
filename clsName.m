function s = clsName(c)
% CLSName  Spike-class index -> short name used by the S3/S4 STA raw-stack helpers.
%   1 -> 'SS'  (simple spike),  2 -> 'CS'  (complex spike)
names = {'SS','CS'};
s = names{c};
end
