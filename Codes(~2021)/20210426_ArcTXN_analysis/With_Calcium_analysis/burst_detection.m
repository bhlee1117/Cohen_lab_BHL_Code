function [Raw IBI_pooled  Pow_spec_pooled Brst_pooled theta_Brst_pooled]=burst_detection(dat,void_th,D_th,post_ref)
max_lag=300;
   if post_ref
    cd_max=2;
   else
    cd_max=4;
   end
   
for cd=1:cd_max
    IBI_pooled{cd}=[];
    Brst_pooled{cd}=[];
    theta_Brst_pooled{cd}=[];
    Pow_spec_pooled{cd}=[];
end
for m=1:size(dat,2)
    Brst_d_pool{m}=[];
    theta_Brst_d_pool{m}=[];
    Pow_spec_d_pooled{m}=[];
    IBI_d_pool{m}=[];
    for d1=1:3
        for c1=1:size(dat{m}.S,1)
            if size(dat{m}.S{c1,d1},2)>1
                moS=movsum(dat{m}.S{c1,3+d1},2+void_th);
                [val loc]=findpeaks(moS);
                S=zeros(1,size(dat{m}.S{c1,3+d1},2)); S(loc)=val;
                S_th=double(S>D_th); S_th(isnan(dat{m}.S{c1,d1}))=NaN;
                Raw.Brst{m,d1}(c1,1)=sum(S_th,'omitnan')/(sum(~isnan(S_th))/30);
                if isempty(find(val>D_th))
                    Raw.Brst{m,d1}(c1,2:3)=NaN;
                else
                    Raw.Brst{m,d1}(c1,2:3)=[mean(val(find(val>D_th))) max(val(find(val>D_th)))];
                end
                total_brst=sum(S_th,'omitnan');
                for lag=1:max_lag
                    Raw.IBI{m,d1}(c1,lag)=sum(S_th(1:end-lag).*S_th(1+lag:end));
                end
                %IBI_prob{m,d1}(c1,:)=IBI{m,d1}(c1,:)/total_brst;
                [acf_tmp l]=autocorr(S_th,max_lag);
                Raw.IBI_prob{m,d1}(c1,:)=acf_tmp(2:end);
%                 [Raw.f Raw.Power_spec{m,d1}(c1,:)]=run_fft(l/30,acf_tmp,0);
                 [F F_tmp]=run_fft([find(~isnan(S_th))-1]/30,S_th(find(~isnan(S_th))),0);
                 Raw.Power_spec{m,d1}(c1,:)=interp1(F,F_tmp(1:size(F,2)),[0:15/6000:15]);    
%                 [Raw.f Raw.Power_spec{m,d1}(c1,:)]=F_tmp;
                Raw.theta_Brst{m,d1}(c1,1)=sum(Raw.IBI{m,d1}(c1,3:5),'omitnan')/(sum(~isnan(S_th))/30);
            else
                Raw.IBI{m,d1}(c1,:)=NaN(1,max_lag);
                Raw.IBI_prob{m,d1}(c1,:)=NaN(1,max_lag);
                Raw.Power_spec{m,d1}(c1,:)=NaN(1,size([0:15/6000:15],2));
                Raw.f=[0:15/6000:15];
                %Raw.Power_spec{m,d1}(c1,:)=NaN(1,max_lag/2+2);
                Raw.Brst{m,d1}(c1,1:3)=NaN(1,3);
                Raw.theta_Brst{m,d1}(c1,1)=NaN;
            end
        end
        IBI_d_pool{m}=[IBI_d_pool{m}; Raw.IBI_prob{m,d1}];
        Brst_d_pool{m}=[Brst_d_pool{m}; Raw.Brst{m,d1}];
        Pow_spec_d_pooled{m}=[Pow_spec_d_pooled{m}; Raw.Power_spec{m,d1}];
        theta_Brst_d_pool{m}=[theta_Brst_d_pool{m};  Raw.theta_Brst{m,d1}];
    end
   if post_ref
    Arc_class{m}=[dat{m}.Arc_post_ref(:,2:end)];
    Arc_class{m}=Arc_class{m}(:);
    cd_max=2;
    cmap=[0.1 0.1 0.1; 1 0 0];
   else
    Arc_class{m}=[dat{m}.Arc(:,2:end)];
    Arc_class{m}=Arc_class{m}(:);
    cd_max=4;
    cmap=distinguishable_colors(4);
   end
    for cd=1:cd_max
        list{m,cd}=Arc_class{m}==cd;
        IBI_pooled{cd}=[IBI_pooled{cd}; IBI_d_pool{m}(find(list{m,cd}),:)];
        Pow_spec_pooled{cd}=[Pow_spec_pooled{cd}; Pow_spec_d_pooled{m}(find(list{m,cd}),:)];
        Brst_pooled{cd}=[Brst_pooled{cd}; Brst_d_pool{m}(find(list{m,cd}),:)];
        theta_Brst_pooled{cd}=[theta_Brst_pooled{cd}; theta_Brst_d_pool{m}(find(list{m,cd}),:)];
    end
    
    % Matrix FR : cell array consist of n X 3 matrix, each column represents
    % days.
    % Quantity of interest is arranged in mouse and groups.
end
%%
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],'Renderer','painters','position',[100 100 400 300]);

for i=1:2
    lineProps.col{1}=cmap(i,:);
    mseb([1:size(IBI_pooled{i},2)]/30,mean(IBI_pooled{i},1,'omitnan'),std(IBI_pooled{i},0,1,'omitnan')./sqrt(sum(~isnan(IBI_pooled{i}),1)),lineProps,1);
    hold all
end
ylabel('Burst probability','FontName','arial rounded mt bold','FontSize',13)
xlabel('Inter-burst interval (s)','FontName','arial rounded mt bold','FontSize',13)
xlim([0 1])
ylim([0 0.1])
legend({'Arc^-','Arc^+'})


figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],'Renderer','painters','position',[100 100 400 300]);

for i=1:2
    lineProps.col{1}=cmap(i,:);
    mseb([0:Raw.f(2):Raw.f(2)*(size(Pow_spec_pooled{i},2)-1)],mean(Pow_spec_pooled{i},1,'omitnan'),std(Pow_spec_pooled{i},0,1,'omitnan')./sqrt(sum(~isnan(Pow_spec_pooled{i}),1)),lineProps,1);
    hold all
end
ylabel('Power (A.U.)','FontName','arial rounded mt bold','FontSize',13)
xlabel('Frequency (Hz)','FontName','arial rounded mt bold','FontSize',13)
xlim([0 15])
ylim([0 0.1])
legend({'Arc^-','Arc^+'})

yaxis_name={'Mean burst rate (Hz)','Mean spike number (/burst)','Max spike number (/event)','Theta burst rate (s^{-1})'};
ylims=[1 10 50];
for i=1:3
    plot_errorbar({Brst_pooled{1}(:,i),Brst_pooled{2}(:,i)},'ranksum',ylims(1,i),yaxis_name{i},{'Arc^-','Arc^+'},0)
end

plot_errorbar({theta_Brst_pooled{1},theta_Brst_pooled{2}},'ranksum',0.1,yaxis_name{4},{'Arc^-','Arc^+'},0)
end