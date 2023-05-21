function P=spike2binnedarray(S)
bin=0.3; %sec;
lag=bin*30;
if size(S,2)>1
tmp=movsum(S,lag,'omitnan');
tmp_isnan=(movsum(isnan(S),lag)>0);
P=double(tmp([lag:lag:length(tmp)])>0);
P(tmp_isnan([lag:lag:length(tmp)]))=NaN;
else
    P=NaN;
end
end