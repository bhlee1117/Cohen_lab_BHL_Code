function perispike_ind=get_perispikeIndex(spike,nTau)

% generate a index matrix, 2024.11.15 Byung Hun Lee
% example: spike = [0 0 0 1 0 0 0]; nTau=[-2:2]
% perispike_ind = [2 3 4 5 6];

perispike_ind=(find(spike==1)'+nTau);
spike_out=sum(perispike_ind<=0 | perispike_ind>length(spike),2)>0;
perispike_ind(spike_out,:)=[];

end