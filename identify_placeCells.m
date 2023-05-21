%function [place_fr field is_pc running]=identify_placeCells(traces,spike,bin,Ard_data,peak_th,in_out_th)
winG=gausswin(5,1); winG=winG/sum(winG);

mean_S=filter(winG,1,traces2place(spike,bin,Ard_data));
mean_F=filter(winG,1,traces2place(traces,bin,Ard_data));
bin_mm=max(Ard_data(:,2))/bin; %mm




    end
end


