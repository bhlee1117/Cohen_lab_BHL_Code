function Sumpeak=extrack_sigma_peak(Peak)
if isempty(Peak)
    Sumpeak=0;
elseif isnan(Peak)
    Sumpeak=NaN; else
    Sumpeak=sum(Peak(:,2)); end 
end