function dat=identiy_pc_batch(dat,Full_result)
day_end=3;
for m=1:size(dat,2)
    for i=1:size(dat{m}.Cal,1)
        for d=1:3
            clear P
            if size(dat{m}.S{i,d},2)>1
                P(:,3)=find(dat{m}.S{i,d}>0)'/30;
                [dat{m}.place_FM_dec{i,d} dat{m}.place_field_dec{i,d} dat{m}.isPC_dec{i,d} R]=identify_place_cell(dat{m}.Cal{i,day_end+d},20,P,Full_result{m}.VR{d},0.12,2.2);
            else
                dat{m}.place_FM_dec{i,d}=NaN;
                dat{m}.place_field_dec{i,d}=NaN;
                dat{m}.isPC_dec{i,d}=NaN;
            end
        end
    end
end
end