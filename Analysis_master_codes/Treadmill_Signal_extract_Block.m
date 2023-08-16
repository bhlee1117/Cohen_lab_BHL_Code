clear
[fpath] = uigetfile_n_dir;
save_name='pcResult_20230716.mat';
%%
time_size=150000;
block_size=15;
DAQ_rate=0.000005;

for i=1:length(fpath)
    load(fullfile(fpath{i},'mcTrace.mat'))
    Result{i}.centers=CellCoordinate;
    Result{i}.fpath=fpath{i};
    disp(['loading  ' fpath{i}])
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[3000:4000]));
Result{i}.FOV=mean(mov_test,3);

    try
        Result{i}.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    end
    try
        Result{i}.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
    end

    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    Result{i}.frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    t_seg=[[1:time_size:length(CamTrigger)] length(CamTrigger)+1];
    t_seg=[t_seg(1:end-1)' t_seg(2:end)'+100];
    t_seg(end)=length(CamTrigger);

    Result{i}.traces=[];
    Result{i}.traces_res=[];
    Result{i}.im_corr=[];
    Result{i}.mcTrace=[];

    for n=1:size(Result{i}.centers,1)


        mov=double(readBinMov([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(1,'%02d') '.bin'], ...
            bs(n)*2+1,bs(n)*2+1));
        mov_res= mov-mean(mov,3);
        bkg = zeros(2, size(mov,3));
        bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
        bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
        mov_res=SeeResiduals(mov_res,mcTrace_block{n,1});
        mov_res= SeeResiduals(mov_res,bkg,1);


        Result{i}.ref_im(1:2*bs(n)+1,1:2*bs(n)+1,n)=mean(mov,3);
        ref_im_vec=tovec(Result{i}.ref_im(1:2*bs(n)+1,1:2*bs(n)+1,n));
        ref_im_vec=(ref_im_vec-mean(ref_im_vec,1))./std(ref_im_vec,0,1);
        mov_mc_vec=tovec(mov);
        mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
        Result{i}.im_corr{n}=[];
        Result{i}.im_corr{n}=[Result{i}.im_corr{n} sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];


        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=mask_footprint([bs(n)+0.5 bs(n)+0.5],movmean(mov_res(:,:,1000:end-1000),10,3),[],8);
        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=imgaussfilt(Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n),0.6);

        for t=2:size(t_seg,1)
            Result{i}.traces(n,t_seg(t-1,1):t_seg(t-1,2)-101)=-(tovec(mov_res(:,:,1:time_size))'*tovec(Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)))';
            Result{i}.mcTrace(:,t_seg(t-1,1):t_seg(t-1,2)-101,n)=mcTrace_block{n,t-1}(:,1:time_size);

            mov=double(readBinMov([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(t,'%02d') '.bin'], ...
                bs(n)*2+1,bs(n)*2+1));

            mov_res= mov-mean(mov,3);
            bkg = zeros(2, size(mov,3));
            bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
            bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
            mov_res=SeeResiduals(mov_res,mcTrace_block{n,t});
            mov_res=SeeResiduals(mov_res,mcTrace_block{n,t}.^2);
            mov_res= SeeResiduals(mov_res,bkg,1);

            mov_mc_vec=tovec(mov);
            mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
            Result{i}.im_corr{n}=[Result{i}.im_corr{n} sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];

        end
        Result{i}.traces(n,t_seg(end,1):t_seg(end,2))=-(tovec(mov_res)'*tovec(Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)))';
        Result{i}.mcTrace(:,t_seg(end,1):t_seg(end,2),n)=mcTrace_block{n,end};


        mcT=squeeze(Result{i}.mcTrace(:,:,n));
        Result{i}.traces_res(n,:)=squeeze(SeeResiduals(reshape(Result{i}.traces(n,:),1,1,[]),mcT));
        Result{i}.traces_res(n,:)=squeeze(SeeResiduals(reshape(Result{i}.traces_res(n,:),1,1,[]),mcT.^2'));
        Result{i}.traces_res(n,:)=squeeze(SeeResiduals(reshape(Result{i}.traces_res(n,:),1,1,[]),mcT(1,:).*mcT(2,:)));
        Result{i}.AvgImg=ave_im{}
    end

end

for i=1:length(fpath)
    Result{i}.spike=zeros(size(Result{i}.traces));
    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),250,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,5.5,2))>0;
end

save(save_name,'Result','fpath','-v7.3')

%%

lap_dist=115;
figure;
ax1=[];
for i=1:length(Result)

    fileList = dir(fullfile(fpath{i}, '*.data'));
    if ~isempty(fileList)
        

        if length(fileList)==1
            fid = fopen(fullfile(fpath{i},fileList.name));
            VRdata = fread(fid,[12 inf],'double');
        else
            error('There is more than one .data file');
        end

        load([fpath{i} '/output_data.mat'])

        %take VR time

        VRdata=VRdata(:,find(VRdata(10,:)>0));
        t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
        t_VR= t_VR-t_VR(1);
        t_VR= milliseconds(t_VR)/1000;
        CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
        CamTrig=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
        Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
        Reward=rescale(Reward(CamTrig))>0.5;

        CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
        t_DAQ=[CamTrig/CamDAQ_rate];
        t_DAQ_scaled=rescale(t_DAQ)*t_VR(end);

        % Calculate converting bin time
        [t_VR arg]=unique(t_VR);
        VRdata=VRdata(:,arg);
        %spike=spike(arg);
        rate=(t_VR(end)-t_VR(1))/(length(t_DAQ_scaled)-1);
        disp(['Calculated rate is ',num2str(rate) ,'sec/frame'])
        new_bin=[t_VR(1):rate:t_VR(end)];

        % Lap change
        lap_end=[0; find(abs(VRdata(8,2:end)-VRdata(8,1:end-1))>0)'; size(VRdata,2)];
        laps=[lap_end(1:end-1)+1 lap_end(2:end)];

        %Cumulative track
        cumTrack=[];
        VRdata(5,:)=VRdata(5,:)-min(VRdata(5,:));
        cumTrack=[cumTrack VRdata(5,laps(1,1):laps(1,2))];

        for l2=2:size(laps,1)
            cumTrack=[cumTrack VRdata(5,laps(l2,1):laps(l2,2))+(l2-1)*lap_dist];
        end
        VRdata(end+1,:)=cumTrack;

        % interpolate
        Virmen_data_int(1,:)=new_bin;
        for j=2:size(VRdata,1)
            Virmen_data_int(j,:)=interp1(t_VR,VRdata(j,:),new_bin,'linear');
        end
        Virmen_data_int(5,:)=Virmen_data_int(end,:);

        ll=1; laps_int=[1];
        lap_trace=zeros(1,size(Virmen_data_int,2));
        % calculate track back from cumulative track
        while ~isempty(find(Virmen_data_int(5,:)>lap_dist))
            sub_ind=find(Virmen_data_int(5,:)>lap_dist);
            laps_int(ll,2)=sub_ind(1);
            lap_trace(laps_int(ll,1):laps_int(ll,2))=ll;
            Virmen_data_int(5,sub_ind)=Virmen_data_int(5,sub_ind)-lap_dist;
            ll=ll+1;
            laps_int(ll,1)=sub_ind(1)+1;
        end
        laps_int(ll,2)=size(Virmen_data_int,2);
        lap_trace(laps_int(ll-1,2)+1:end)=ll;

        Virmen_data_int(8,:)=round(Virmen_data_int(8,:));
        Result{i}.Virmen=Virmen_data_int;
        nexttile([1 1])
        plot(Result{i}.Virmen(1,:),rescale(Result{i}.Virmen(5,:)),Result{i}.Virmen(1,:),Reward)
        hold all
        line(Result{i}.Virmen(1,[1 end]),[0.8 0.8],'color','r')
        legend('Virmen Track','Reward at DAQ')
        title(fpath{i},'Interpreter','none')

    end

end
save(save_name,'Result','fpath','-v7.3')
