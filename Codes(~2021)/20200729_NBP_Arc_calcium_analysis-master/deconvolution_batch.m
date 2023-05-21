function dat_dc_added=deconvolution_batch(dat)
dat_dc_added=dat;
for m=1:size(dat,2)
    m
    for i=1:size(dat{m}.Cal,1)
        i
        for d=1:3
            if size(dat{m}.Cal{i,d},2)==1
                dat_dc_added{m}.S{i,d}=NaN;
                dat_dc_added{m}.S{i,d+3}=NaN;
            else
                S=deconvolve_BH(dat{m}.Cal{i,d},20,200,0.1); % ms
                dat_dc_added{m}.S{i,d}=S(1:size(dat{m}.run{d},1));
                dat_dc_added{m}.S{i,d+3}=dat_dc_added{m}.S{i,d};
                dat_dc_added{m}.S{i,d+3}(isnan(dat{m}.Cal{i,3+d}(1:size(dat{m}.run{d},1))) | (dat{m}.run{d}==0)')=NaN;
                dat_dc_added{m}.S{i,d+6}=dat_dc_added{m}.S{i,d};
                dat_dc_added{m}.S{i,d+6}(isnan(dat{m}.Cal{i,3+d}(1:size(dat{m}.run{d},1))))=NaN;
                
            end
        end
    end
end
end