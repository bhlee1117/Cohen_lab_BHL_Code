function alp_patterns = dragonflyHadamardPatterns_yq(dmd,m_q,dmd_pat) % hadamard_mode)
%
%   2018-2020 Vicente Parot/Yitong Qi
%   Cohen Lab - Harvard university


nstims = 1;
blocksize = m_q;
elementsize = 1;
alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize); %theoretical pattern

% prepare for actual acquisition
alp_patterns_cam = cat(3,alp_patterns,alp_patterns(:,:,end)); 
alp_patterns_lt = eye(size(alp_patterns_cam,3)*2);
alp_patterns_lt(1:2:end,:) = [];
alp_patterns_cam_all =  double(reshape(alp_patterns_cam,numel(alp_patterns_cam(:,:,1)),size(alp_patterns_cam,3))) *...
                        alp_patterns_lt;
alp_patterns_cam_all = uint8(reshape(alp_patterns_cam_all,size(alp_patterns_cam,1),size(alp_patterns_cam,2),[]));

figure;moviesc(alp_btd_to_logical(alp_patterns_cam_all))

alp_patterns = alp_logical_to_btd(permute(bsxfun(@and,alp_btd_to_logical(repmat(alp_patterns_cam_all,[1 1 nstims])), dmd_pat),[2 1 3]));
% alp_patterns = alp_logical_to_btd(bsxfun(@and,alp_btd_to_logical(repmat(alp_patterns_cam_all,[1 1 nstims])), dmd_pat));



end
