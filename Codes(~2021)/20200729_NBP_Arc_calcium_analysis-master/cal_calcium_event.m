function [datum cond_Peaks_st result]=cal_calcium_event(Full_result,mouse,day,cond_day)
comb=[1 1;1 2;2 1;2 2];
cmap=[0.5 0.5 0.5 ;distinguishable_colors(3)];
%comb_day=[1 1 1;2 1 1;1 2 1;1 1 2;2 2 1;2 1 2;1 2 2;2 2 2];
comb_day=[2];
for m=mouse
    
    
    for i=day{m}  % Day of interest (may different mouse by mouse)
        datum{m}.Arc=Full_result{m}.Arc_class(Full_result{1, 1}.list_identified(:,1),4:end);
        
        for j=1:size(Full_result{1, 1}.list_identified,1)
            for c=1:size(comb,1)
                if datum{m}.Arc(j,2*i-1)==comb(c,1) && datum{m}.Arc(j,2*i)==comb(c,2)
        datum{m}.Arc_class(j,i)=c;
                end
            end
            if isnan(Full_result{1, 1}.list_identified(j,i+1))  % Not detected
                
        datum{m}.cal_sigma{j,i}=[]; % Sigma
        datum{m}.dFF{j,i}=[]; % dF/F
        datum{m}.Peaks{j,i}=nan(1,5);
            else
        datum{m}.cal_sigma{j,i}=full(Full_result{1, 1}.Calcium{1,i}.cal_sigma(Full_result{1, 1}.list_identified(j,i+1),:));
        datum{m}.dFF{j,i}=full(Full_result{1, 1}.Calcium{1,i}.C_df(Full_result{1, 1}.list_identified(j,i+1),:));
        datum{m}.time{j,i}=size(datum{m}.dFF{j,i},2)*1/30;
        datum{m}.Peaks{j,i}=peak_finding(datum{m}.cal_sigma{j,i},datum{m}.dFF{j,i},...
                                    Full_result{1, 1}.Calcium{1,i}.ini_fin{1,Full_result{1, 1}.list_identified(j,i+1)});  
            end
        end
    end
end

for m=mouse
        for c=1:size(comb,1)
            [rc cc]=find(datum{m}.Arc_class==c);
            for g=1:size(rc,1)
                if isempty(datum{m}.Peaks{rc(g,1),cc(g,1)})
                ave_SP{c}(g,:)=zeros(1,5);    
                else
                    if isnan(datum{m}.Peaks{rc(g,1),cc(g,1)}(1,1))
                ave_SP{c}(g,:)=nan(1,5);        
                    else
                ave_SP{c}(g,:)=[sum(datum{m}.Peaks{rc(g,1),cc(g,1)}(:,1:2),1) ...
                                mean(datum{m}.Peaks{rc(g,1),cc(g,1)}(:,1:2),1) ...
                                size(datum{m}.Peaks{rc(g,1),cc(g,1)},1)./datum{m}.time{rc(g,1),cc(g,1)}];
                    end
                end
            end
        end
        
        for cd=1:size(comb_day,1)
            g=1;
            for d=cond_day{m}
            l(:,g)=datum{m}.Arc_class(:,d)==comb_day(cd,d);
            g=g+1;
            end
            for j=find(min(l,[],2))'
            for d=day{m}
            cond_Peaks{cd}{j,d}=datum{m}.Peaks{j,d};
            if isempty(datum{m}.Peaks{j,d})        % No peak
            cond_Peaks_st{cd}{j,d}=zeros(1,5);
            else                               
                if isnan(datum{m}.Peaks{j,d}(1,1)) % Spatial component not detected
            cond_Peaks_st{cd}{j,d}=nan(1,5);        
                else
            cond_Peaks_st{cd}{j,d}=[sum(datum{m}.Peaks{j,d}(:,1:2),1) mean(datum{m}.Peaks{j,d}(:,1:2),1)...
                                   size(datum{m}.Peaks{j,d},1)./datum{m}.time{j,d}];
                end
            end
            g=g+1;
            end
            end
            [s1 s2]=cellfun(@size,cond_Peaks_st{cd});
            g2=1;
            for kk=find(s1(:,1))'
            result{m,cd}(g2,:)=cell2mat(cond_Peaks_st{cd}(kk,:));
            g2=g2+1;
            end
        end
        
        
end
% 
% 
% 
% 
% end