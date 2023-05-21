%load('\\Neurobiophysics\Byunghun_Lee\Reports\2021\20210330_RevisionÀÚ·á\Figure_simulation\20210410_response_ftn.mat')
time_bin=1; %min
period=60*50; %min
CFC_time= 60*24;
basal_level=0.5;
%reponse_ftn_TXN=[0  11.4606   21.3767   29.6686   36.2566   41.0609   44.0019   45.0000   44.4231   43.0000   38.6169   32.0000   24.8296   18.0000   12.6929    9.0000    7.4375    6.0000    3.4375    0];
%reponse_ftn=int_n;
% response_ftn_TXN=intp_dat;
% response_ftn=x2;
%%
N=500;
%basal=[0.5:-0.02:0];
%% Variation in threshold
threshold=0.2;
clear Conv_trac   A_signal_trace sum_A_trace sum_conv_trace
% b=1;
threshold=[0.15:0.2:0.95];
for b=1:5
    for it=1:10
        [C A]=activitiy_generator(N,1/10);
        Conv_trac(it,2*b-1:2*b)=C; A_signal_trace{it,b}=A;
        %     Conv_trac{it,2*b-1}(Conv_trac{it,2*b-1}>1)=1;
        %     Conv_trac{it,2*b}(Conv_trac{it,2*b}>1)=1;
        Conv_trac_bin{it,2*b-1}=Conv_trac{it,2*b-1}>threshold(b);
        Conv_trac_bin{it,2*b}=Conv_trac{it,2*b}>threshold(b);
        %
%         TP{b}(it,1)=sum(Conv_trac{it,2*b-1}(:,26*60)>0.9 & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
%         %sum(Conv_trac{it,2*b-1}(:,26*60)>0.9);
%         TP{b}(it,2)=sum(max(Conv_trac{it,2*b}(:,24*60+3:24*60+7),[],2)>0.2 & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
%         %sum(max(Conv_trac{it,2*b}(:,24*60+5),[],2)>0.2);
%         TN{b}(it,1)=sum(Conv_trac{it,2*b-1}(:,26*60)<0.9 & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
%         TN{b}(it,2)=sum(max(Conv_trac{it,2*b}(:,24*60+3:24*60+7),[],2)<0.2 & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
%         FN{b}(it,1)=sum(Conv_trac{it,2*b-1}(:,26*60)<0.9 & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
%         FN{b}(it,2)=sum(max(Conv_trac{it,2*b}(:,24*60+3:24*60+7),[],2)<0.2 & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
%         FP{b}(it,1)=sum(Conv_trac{it,2*b-1}(:,26*60)>0.9 & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
%         FP{b}(it,2)=sum(max(Conv_trac{it,2*b}(:,24*60+3:24*60+7),[],2)>0.2 & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
%         Ac{b}(it,1)=(sum(Conv_trac{it,2*b-1}(:,26*60)>0.9 & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2))...
%             + sum(Conv_trac{it,2*b-1}(:,26*60)<0.9 & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2)))/N;%/...
%         
%         Ac{b}(it,2)=(sum(Conv_trac{it,2*b}(:,24*60+5)>0.2 & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2))...
%             +sum(Conv_trac{it,2*b}(:,24*60+5)<0.2 & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2)))/N;%/...
       
    end
    sum_conv_trace{b}=cellfun(@sum,Conv_trac_bin(:,2*b-1:2*b),'UniformOutput',false);
    sum_A_trace{b}=cell2mat(cellfun(@sum,A_signal_trace(:,b)','UniformOutput',false)')/N;
end
%% plot
Wang =[0 30 60 120 240 720;2 34 63 80 51 8];
    figure
    cmap=distinguishable_colors(5);
for b=1:5
    M_GFP=mean(cell2mat(sum_conv_trace{b}(:,1)))/N; S_GFP=std(cell2mat(sum_conv_trace{b}(:,1)))/N;
    lineProps.col{1}=cmap(b,:);
    mseb([1:1:size(M_GFP,2)]/60,M_GFP,S_GFP,lineProps,1);
    hold all
end
    set(gca,'Xtick',[22 24 28 32 36],'XtickLabel',[-2 0 4 8 12])
    xlabel('Time (hour)')
    ylabel('Fraction of Arc+ neuron')
    xlim([22 36])
    ylim([0 0.6])
figure
for b=1:5
    M_TXN=mean(cell2mat(sum_conv_trace{b}(:,2)))/N;  S_TXN=std(cell2mat(sum_conv_trace{b}(:,2)))/N;
    lineProps.col{1}=cmap(b,:);
    mseb([1:1:size(M_GFP,2)]/60,M_TXN,S_TXN,lineProps,1);
    hold all
end
    xlabel('Time (min)')
    ylabel('Fraction of Arc+ neuron')
    set(gca,'Xtick',[23.75 24 24.25 24.5],'XtickLabel',[-15 0 15 30])
    xlim([23.7 24.5])
    ylim([0 0.3])
    %% Variation in basal level
basal=[1]/10;
th=[0.95 0.15];
clear Conv_trac   A_signal_trace sum_A_trace sum_conv_trace TP FP FN TN
for b=1:length(basal)
    for it=1:10
        [C A]=activitiy_generator(N,basal(b));
        Conv_trac(it,2*b-1:2*b)=C; A_signal_trace{it,b}=A;
        %     Conv_trac{it,2*b-1}(Conv_trac{it,2*b-1}>1)=1;
        %     Conv_trac{it,2*b}(Conv_trac{it,2*b}>1)=1;
        Conv_trac_bin{it,2*b-1}=Conv_trac{it,2*b-1}>th(1);
        Conv_trac_bin{it,2*b}=Conv_trac{it,2*b}>th(2);
        %
        TP{b}(it,1)=sum(max(Conv_trac{it,2*b-1}(:,25.5*60:26.5*60),[],2)>th(1) & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
        TP{b}(it,2)=sum(max(Conv_trac{it,2*b}(:,24*60+4:24*60+7),[],2)>th(2) & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
        
        TN{b}(it,1)=sum(max(Conv_trac{it,2*b-1}(:,25.5*60:26.5*60),[],2)<th(1) & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
        TN{b}(it,2)=sum(max(Conv_trac{it,2*b}(:,24*60+4:24*60+7),[],2)<th(2) & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
        
        FN{b}(it,1)=sum(max(Conv_trac{it,2*b-1}(:,25.5*60:26.5*60),[],2)<th(1) & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
        FN{b}(it,2)=sum(max(Conv_trac{it,2*b}(:,24*60+4:24*60+7),[],2)<th(2) & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
        
        FP{b}(it,1)=sum(max(Conv_trac{it,2*b-1}(:,25.5*60:26.5*60),[],2)>th(1) & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
        FP{b}(it,2)=sum(max(Conv_trac{it,2*b}(:,24*60+4:24*60+7),[],2)>th(2) & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2));
        
        Ac{b}(it,1)=(sum(max(Conv_trac{it,2*b-1}(:,25.5*60:26.5*60),[],2)>th(1) & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2))...
            + sum(max(Conv_trac{it,2*b-1}(:,25.5*60:26.5*60),[],2)<th(1) & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2)))/N;%/...
        Ac{b}(it,2)=(sum(max(Conv_trac{it,2*b}(:,24*60+4:24*60+7),[],2)>th(2) & max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2))...
                        +sum(max(Conv_trac{it,2*b}(:,24*60+4:24*60+7),[],2)<th(2) & ~max(A_signal_trace{it,b}(:,24*60-3:24*60),[],2)))/N;%/...
       
    end
    sum_conv_trace{b}=cellfun(@sum,Conv_trac_bin(:,2*b-1:2*b),'UniformOutput',false);
    sum_A_trace{b}=cell2mat(cellfun(@sum,A_signal_trace(:,b)','UniformOutput',false)')/N;
end



%% Calculate errors
a=cellfun(@mean,FP,'UniformOutput',false);
a=cell2mat(a')/N
b=cellfun(@std,FP,'UniformOutput',false);
b=cell2mat(b')+1/N
%%
figure
lineProps.col{1}=[0 0.2 1];
mseb(flipud(basal)*10,a(:,2)',b(:,2)',lineProps,1);
lineProps.col{1}=[1 0.2 0];
hold all
mseb(flipud(basal)*10,a(:,1)',b(:,1)',lineProps,1);
xlabel('Basal level (% per 10 min)')
ylabel('Accuracy')
%set(gca, 'XDir','reverse')
