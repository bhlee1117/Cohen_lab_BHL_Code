function dat_f=fourier_batch(dat)
dat_f=dat;
for m=1:size(dat,2)
    for c1=1:size(dat{m}.S,1)
        for d1=1:3
            if size(dat{m}.S{c1,d1},2)>1
             [time_vec FT]=run_fft([find(~isnan(dat{m}.S{c1,3+d1}))-1]/30,...
                             dat{m}.S{c1,3+d1}(find(~isnan(dat{m}.S{c1,3+d1}))),0);
                    dat_f{m}.FT{c1,3+d1}=interp1(time_vec,FT(1:size(time_vec,2)),[0:15/6000:15]);
                    
              [time_vec FT]=run_fft([find(~isnan(dat{m}.S{c1,d1}))-1]/30,...
                             dat{m}.S{c1,d1}(find(~isnan(dat{m}.S{c1,d1}))),1);
                    dat_f{m}.FT{c1,d1}=interp1(time_vec,FT(1:size(time_vec,2)),[0:15/6000:15]);    
            else
                  dat_f{m}.FT{c1,d1}=NaN; dat_f{m}.FT{c1,3+d1}=NaN;
            end
        end
    end
end
end