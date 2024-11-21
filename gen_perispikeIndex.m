function perispike_ind=gen_perispikeIndex(spike,nTau)

% generate a index matrix,
% example: spike = [0 0 0 1 0 0 0]; nTau=[-2:2]
% perispike_ind = [2 3 4 5 6];

perispike_ind=(find(spike)'+nTau);
spike_out=sum(perispike_ind<=0 | perispike_ind>length(spike));
perispike_ind(spike_out,:)=[];

end