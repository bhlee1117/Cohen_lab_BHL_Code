function fr=decon_sum(S)
fr=sum(S,'omitnan')/sum(~isnan(S))*30;
end