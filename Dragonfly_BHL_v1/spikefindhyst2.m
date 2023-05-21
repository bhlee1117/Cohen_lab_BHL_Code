function out = spikefindhyst2(dat, thresh1, thresh2, min_spk_width);
% Function out = spikefindhyst(dat, thresh1, thresh2);
% for each time that dat goes above thresh1 and then below thresh2, returns the index of the local
% maximum in dat.  Should have thresh1 > thresh2.
% Use of different thresholds on the rising and falling edges guards
% against spurious triggers from noise near threshold.
% DRH and AEC 24 Feb. 2011
% AEC 10 May 2015

downstroke = 1;
out = [];
while 1
    upstroke = find(dat(downstroke:end) > thresh1, 1, 'first');
    if isempty(upstroke)
        break
    else;
        upstroke = upstroke + downstroke - 1;
    end;
    
    downstroke = find(dat(upstroke:end) < thresh2, 1, 'first');
    
    if isempty(downstroke)
        break
    else 
        downstroke = downstroke + upstroke - 1;
    end;
    
    if (downstroke-upstroke)<min_spk_width
        continue
    end
    [~,indx] = max(dat(upstroke:downstroke));
    out = [out (indx + upstroke - 1 )];
end;

return