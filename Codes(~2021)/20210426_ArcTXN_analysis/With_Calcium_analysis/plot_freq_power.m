function plot_freq_power(HR,rng,rep_d,method)

    for d1=rep_d
    dd_tmp=HR(:,d1)';
    for g=1:size(dd_tmp,2)
    dddd{g}=sum(dd_tmp{g}(:,rng),2);
    end
    plot_errorbar(dddd,method,1,'Prob.',{'Arc^-','Arc^+'},1)
    end
end