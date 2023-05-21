%% This code is written for analyzing Calcium imaging + Arc-TXN population 
% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy,
%           Seoul National University, 2020/09/17
%           Add place field identification, 2020/09/20
%           Add classify_running funciton with very generous threshold,
%           2020/12/29
function dat=show_cdf_arc(Full_result,mouse)
clear dat
comb=[1 1;1 2;2 1;2 2];
for m=mouse
    day_end=(size(Full_result{m}.Arc_class,2)-3)/2;
    dat{m}.Arc(:,1)=Full_result{m}.list_identified(:,1);
    
        for i=1:day_end % Day of interest
         for k=1:size(comb,1) % Classification in Arc TXN
             l=find(Full_result{m}.Arc_class(dat{m}.Arc(:,1),2*i+2)==comb(k,1) &...
                    Full_result{m}.Arc_class(dat{m}.Arc(:,1),2*i+3)==comb(k,2));    % find the list (identified list standard) of each combination
             dat{m}.Arc(l,i+1)=k;
         end
        dat{m}.Arc(dat{m}.Arc==0)=5; % Na 
        dat{m}.Arc_post_ref=dat{m}.Arc;
        dat{m}.Arc_post_ref(dat{m}.Arc==1 | dat{m}.Arc==3)=1;
        dat{m}.Arc_post_ref(dat{m}.Arc==2 | dat{m}.Arc==4)=2;
        end
                    m
        for i=1:size(Full_result{m}.list_identified,1) % Cells in the identified list
            i
            for d=1:day_end
        cell=Full_result{m}.list_identified(i,d+1); % i th (Arc reference) -> Calcium matched cell # in day
        if isnan(cell)
        dat{m}.Cal{i,d}=NaN;
        dat{m}.Cal{i,3+d}=NaN;    
        dat{m}.Cal{i,6+d}=NaN;    
        dat{m}.Cal{i,9+d}=NaN;   
        dat{m}.Peak{i,d}=NaN;
        dat{m}.place_FM{i,d}=NaN;
        dat{m}.place_field{i,d}=NaN; 
        dat{m}.isPC{i,d}=NaN;
        dat{m}.Tr{i,d}=NaN;

        else
        dat{m}.Cal{i,d}=full(Full_result{m}.Calcium{d}.C_df(cell,:));              % C_df
        dat{m}.Cal{i,day_end+d}=full(Full_result{m}.Calcium{d}.C_df(cell,:));      % off focus -> NaN
        dat{m}.Cal{i,2*day_end+d}=full(Full_result{m}.Calcium{d}.ini_fin{1,cell}); % Transients 
        dat{m}.Cal{i,3*day_end+d}=full(Full_result{m}.Calcium{d}.off_focus{1,cell});   
        for foc_out=1:size(dat{m}.Cal{i,3*day_end+d},1)
            dat{m}.Cal{i,day_end+d}(dat{m}.Cal{i,3*day_end+d}(foc_out,1):dat{m}.Cal{i,3*day_end+d}(foc_out,2)-1)=NaN;
        end
        dat{m}.cal_sigma{i,d}=Full_result{m}.Calcium{d}.cal_sigma(cell,:);

        [dat{m}.Peak{i,d} dat{m}.Tr{i,d}]=peak_finding(dat{m}.cal_sigma{i,d},dat{m}.Cal{i,d},dat{m}.Cal{i,2*day_end+d}); % Peaks
        dat{m}.run{d}=classify_running(Full_result{m}.VR{d});
        dat{m}.Cal{i,day_end+d}(1,~dat{m}.run{d})=NaN;
     
        end
            end
        end
end

end