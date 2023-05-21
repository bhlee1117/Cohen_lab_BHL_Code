function alp_patterns = dragonflyRCaMPtopatchPatterns(dmd) % hadamard_mode)
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%
%             roirows = 19*10+0:57*10+3; % extended central quad 2017-03-17
%             roicols = 32*10+0:69*10+5;
%             hadamard_mode = 'voltage';
%             device = [];
% ncols = 1024;
% nrows = 768;
ncols = dmd.device.height;
nrows = dmd.device.width;

super_mask = ones(1,nrows,ncols);
superstim_mask = super_mask;
stim0_mask = ones(nrows,ncols);

% switch hadamard_mode
%     case 'voltage'
%         sequence_mode = 'voltage_stim'; % structural
%     case 'structural'
%         sequence_mode = 'shorter_63c14'; % structural
%     case 'activity'
%         sequence_mode = 'hadamard_11c3_200s_nostim'; % activity
        sequence_mode = 'hadamard_11c3_12periods_pulse'; % had sync
        
            roirows = 18*10+0:55*10+3; % extended central quad 2017-12-18
            roicols = 32*10+0:69*10+5;
% end

% % stim0_mask = zeros(nrows,ncols);
% sequence_mode = 'hadamard_11c3_sprinkledstim'; % had sync
% roirows = 19*10+0:57*10+3; % extended central quad 2017-03-17
% roicols = 32*10+0:69*10+5;

% stimroirows = 33*10+0:57*10+3; % LTP stim roi 2017-04-07
% stimroicols = 32*10+0:69*10+5;


if true
    % full ROI
    super_mask(1,roirows,roicols) = 1;
%     figure; imshow2(squeeze(super_mask));
%     superstim_mask(1,:,:) = 1;
    superstim_mask(1,roirows,roicols) = 1;
else
    % quadrants 
    roirowslime = roirows(1:ceil(end/2));
    super_mask(1,roirowslime,roicols) = 1;
    roicolsblue = roicols(1:ceil(end/2));
    superstim_mask(1,roirows,roicolsblue) = 1;
end
%% stim params

% stim_mode = 'checkerboards';
% period = 1000/9.6; % dmd pixels
% checkerboard = (double(mod(0:nrows-1,period)<period/2)'-.5)*(double(mod(0:ncols-1,period)<period/2)-.5)>0;
% 
% stim_mode = 'bgspots';
bglevel = 0; % 2^-4; % 0; % % 2^-6;
% load('D:\Code\Hadamard runtime\Analysis\201838_stimulation_masks.mat')      % cortex
% cort_mask = stimulation_masks(:,:,1);
% load('D:\Code\Hadamard runtime\Analysis\195716_stimulation_masks.mat')      % cortex
% thal_mask = stimulation_masks(:,:,1);


% stim0_mask = cort_mask + thal_mask;

% stim0_mask = stimulation_masks(:,:,1); % comment for whole field
% stim0_mask = stimpat';

%% pulses
npulses = 8; 
pulselength = 5; % pulses
% pulselength = 50; % blocks
pulseperiod = 50;
totallength = 400;
x = (0:.8:totallength);
spotmod = mod(x,pulseperiod)<pulselength & (totallength-pulseperiod*npulses)<=x & x<totallength;  
% spotmod = circshift(spotmod,300,2);
% figure; plot(x,spotmod) % pulses
%%
%         x = (0:.8:400); spotmod = 5<=x&x<385; % plot(x,spotmod) % dc

super_mask = uint8(super_mask);
stim0_mask = permute(uint8(superstim_mask),[2 3 1]).*uint8(stim0_mask);
selective_mask = ones(nrows,ncols,1);
% selective_mask = new_peaks_mask;
super_mask = super_mask.*permute(uint8(selective_mask),[3 1 2]);
switch sequence_mode
    case 'hadamard_11c3_12periods_pulse'
        ndatapoints = 3; % 11
%         clc
        blocksize = [43 12]; % [11 3]
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,permute(super_mask,[3 2 1]));
%         stimmask = permute(super_mask,[2 3 1]);
        npatsblue = 501;
        stim_spot_only = false(nrows,ncols,npatsblue);
        stim_spot_only(:,:,:) = bsxfun(@or,stim0_mask(:,:),stim_spot_only(:,:,:));
        stim_spot_only(:,:,:) = bsxfun(@and,permute(spotmod,[1 3 2]),stim_spot_only(:,:,:));
        stim_spot_only(:,:,end) = false;
%         close all
%         figure windowstyle docked
%         imshow([mean(stim_spot_only,3)],[])
%         title({'randonly randspot'})
        stim_spot_only = alp_logical_to_btd(stim_spot_only);
        alp_patterns = repmat(alp_patterns,[1 1 ndatapoints]);
        alp_patterns = cat(3,alp_patterns(:,:,1:2)*0,alp_patterns,alp_patterns(:,:,1:2)*0);
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
        alp_patterns = cat(3,alp_patterns,stim_spot_only,alp_patterns);
%         alp_patterns(:,:,17) = 0;
    case 'voltage_stim' % the one that gets executed
        datadir = 'R:';

        sdir = dir(fullfile(datadir,'S*'));
        sdir = sdir(cell2mat({sdir.isdir}));
        [~,ix] = sort(arrayfun(@(x)x.name,sdir,'uni',false));
        sdir = sdir(ix(end));

        ramslicepath = fullfile(datadir,sdir.name);

        fdir = dir(fullfile(ramslicepath,'FOV*'));
        fdir = fdir(cell2mat({fdir.isdir}));
        [~,ix] = sort(arrayfun(@(x)x.name,fdir,'uni',false));
        fdir = fdir(ix(end));

        ramfovpath = fullfile(ramslicepath,fdir.name);
        
        %%
        dl = dir(fullfile(ramfovpath,'*lineclick*'));
        dl = sort({dl.name})';
        assert(~isempty(dl),sprintf('lineclick not found in\n%s',ramfovpath))
        dmd_cal_existing_fname = fullfile(ramfovpath,dl{end});
        load(dmd_cal_existing_fname,'xyzRowColDMD')
%         xyzRowColDMD = xyzRowColDMD - [85 35 0];

        roundRobin = true;
        if roundRobin 
            pat = zeros(ncols,nrows,size(xyzRowColDMD,1)+1);
            for it = 1:size(xyzRowColDMD,1)
                pat(max(min(round(xyzRowColDMD(it,1)),end),1),max(min(round(xyzRowColDMD(it,2)),end),1),it) = 1;
            end
            pat(:,:,end) = min(sum(pat(:,:,1:end-1),3),1);
        else
            pat = zeros(ncols,nrows);
            for it = 1:size(xyzRowColDMD,1)
                pat(max(min(round(xyzRowColDMD(it,1)),end),1),max(min(round(xyzRowColDMD(it,2)),end),1)) = 1;
            end
        end
        pat = imdilate(pat,ones(27)); 
%         pat = circshift(pat,[-90 -40]);
        
%         figure windowstyle docked
%         moviefixsc(flip(rot90(pat)))
        alp_patterns = alp_logical_to_btd(pat);
%         alp_patterns = reshape([alp_patterns alp_patterns],[size(alp_patterns) 1].*[1 1 2]);
    case 'voltage_stim_old'
        datadir = 'R:';

        sdir = dir(fullfile(datadir,'S*'));
        sdir = sdir(cell2mat({sdir.isdir}));
        [~,ix] = sort(arrayfun(@(x)x.name,sdir,'uni',false));
        sdir = sdir(ix(end));

        ramslicepath = fullfile(datadir,sdir.name);

        fdir = dir(fullfile(ramslicepath,'FOV*'));
        fdir = fdir(cell2mat({fdir.isdir}));
        [~,ix] = sort(arrayfun(@(x)x.name,fdir,'uni',false));
        fdir = fdir(ix(end));

        ramfovpath = fullfile(ramslicepath,fdir.name);
        
        dirlist = dir(fullfile(ramslicepath,'*dmd_cal*.mat'));
        [~,ix] = sort(arrayfun(@(x)x.name,dirlist,'uni',false));
%         [~, ix] = sort(cell2mat({dirlist.name}));
%         msgbox(num2str(size(ix)))
        load(fullfile(ramslicepath,dirlist(ix(end)).name));

        dirlist = dir(fullfile(ramfovpath,'*click_spots*.mat'));
%         [~,ix] = sort(arrayfun(@(x)x.name,dirlist,'uni',false));
        [~, ix] = sort(cell2mat({dirlist.datenum}));
%         [~, ix] = sort(cell2mat({dirlist.name}));
        load(fullfile(ramfovpath,dirlist(ix(end)).name));

        dirlist = dir(fullfile(ramfovpath,'*GFP.tif'));
        [~, ix] = sort(cell2mat({dirlist.datenum}));
        click_plimg = imread(fullfile(ramfovpath,dirlist(ix(end)).name));
        [cplr, cplc] = size(click_plimg);
        %%
% [-(lastClickedX+cal_im_c/2-cplc/2) -(lastClickedY+cal_im_r/2-cplr/2) lastClickedX*0+1]
        stimDotsRowColDMD = [lastClickedX+cal_dmd_im_c/2-cplc/2 lastClickedY+cal_dmd_im_r/2-cplr/2 lastClickedX*0+1]*stimCalMatrix;
        stimDotsRowColDMD = round(stimDotsRowColDMD);
        pat = zeros(ncols,nrows);
        for it = 1:size(stimDotsRowColDMD,1)
            pat(max(min(stimDotsRowColDMD(it,1),end),1),max(min(stimDotsRowColDMD(it,2),end),1)) = 1;
        end
        pat = imdilate(pat,ones(281)); 
        
%         pat = zeros(device.height,device.width);
%         pat(roicols,roirows) = 1;

% device.put(pat*255)

        figure windowstyle docked
        imshow(pat)
        %%
        alp_patterns = alp_logical_to_btd(pat);
        alp_patterns = reshape([alp_patterns alp_patterns],[size(alp_patterns) 1].*[1 1 2]);
    case 'hadamard_11c3_spotonly_once'
        ndatapoints = 11;
        clc
        blocksize = [11 3];
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,permute(super_mask,[3 2 1]));
%         stimmask = permute(super_mask,[2 3 1]);
        npatsblue = 501;
        stim_spot_only = false(nrows,ncols,npatsblue);
        stim_spot_only(:,:,:) = bsxfun(@or,stim0_mask(:,:),stim_spot_only(:,:,:));
        stim_spot_only(:,:,:) = bsxfun(@and,permute(spotmod,[1 3 2]),stim_spot_only(:,:,:));
        stim_spot_only(:,:,end) = false;
%         close all
%         figure windowstyle docked
%         imshow([mean(stim_spot_only,3)],[])
%         title({'randonly randspot'})
        stim_spot_only = alp_logical_to_btd(stim_spot_only);
        alp_patterns_pre = cat(3,...
            false(size(alp_patterns(:,:,ones(4,1)))),...
            repmat(alp_patterns,[1 1 ndatapoints]));
        alp_patterns_pre = reshape([alp_patterns_pre alp_patterns_pre*0],size(alp_patterns_pre).*[1 1 2]);
        alp_patterns_post = cat(3,...
            repmat(alp_patterns,[1 1 ndatapoints]),...
            false(size(alp_patterns(:,:,ones(4,1)))));
        alp_patterns_post = reshape([alp_patterns_post alp_patterns_post*0],size(alp_patterns_post).*[1 1 2]);
                alp_patterns = cat(3,alp_patterns_pre,stim_spot_only,alp_patterns_post);
%         a02_load_seq_firefly
%         clear alp_patterns selective_mask stim0_mask stim1_mask stim2_mask stim_bg_nospot stim_spot_only check1 check2 alp_patterns_pre alp_patterns_post
%         disp([sequence_mode  ' loaded'])
    case 'hadamard_11c3_200s_nostim'
        clc
        blocksize = [11 3];
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,permute(super_mask,[3 2 1]));
%         alp_patterns = bsxfun(@times,alp_patterns,permute(uint8([1 0 0 0 0 0 0 0 0 0 0 0]),[1 3 2]));
        close all
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
        figure windowstyle docked
        moviesc(alp_btd_to_logical(alp_patterns))
        disp([sequence_mode  ' loaded'])
    case 'hadamard_11c3_sprinkledstim'
%         nstims = 5;
        ndatapoints = 11;
        clc
        blocksize = [11 3];
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,permute(super_mask,[3 2 1]));
        stimmask = permute(super_mask,[2 3 1]);
        npatsblue = 501;
        stim_bg_nospot = false(nrows,ncols,npatsblue);
        stim_spot_only = false(nrows,ncols,npatsblue);
        stim_bg_nospot(roirows,roicols,:) = rand(numel(roirows),numel(roicols),npatsblue) <= bglevel;
        stim_bg_nospot(roirows,roicols,:) = bsxfun(@and,~stim0_mask(roirows,roicols),stim_bg_nospot(roirows,roicols,:));
        stim_spot_only(roirows,roicols,:) = bsxfun(@or,stim0_mask(roirows,roicols),stim_spot_only(roirows,roicols,:));
        stim_spot_only(roirows,roicols,:) = bsxfun(@and,permute(spotmod,[1 3 2]),stim_spot_only(roirows,roicols,:));
        stim_bg_nospot(:,:,end) = false;
        stim_spot_only(:,:,end) = false;
        close all
        figure windowstyle docked
        imshow([mean(stim_bg_nospot,3) mean(stim_spot_only,3)],[])
        title({'randonly randspot'})
%         check1 = bsxfun(@or, checkerboard & permute(super_mask,[2 3 1]),false(size(stim_bg_nospot)));
%         check1 = bsxfun(@and,permute(spotmod,[1 3 2]),check1);
%         check1 = alp_logical_to_btd(check1);
%         check2 = bsxfun(@or,~checkerboard & permute(super_mask,[2 3 1]),false(size(stim_bg_nospot)));
%         check2 = bsxfun(@and,permute(spotmod,[1 3 2]),check2);
%         check2 = alp_logical_to_btd(check2);
%         check1(:,:,end) = false;
%         check2(:,:,end) = false;
%         check1(:,:,end-12:end) = false;
%         check2(:,:,end-12:end) = false;
        stim_bg_nospot = alp_logical_to_btd(stim_bg_nospot);
        stim_spot_only = alp_logical_to_btd(stim_spot_only);
        alp_patterns_pre = cat(3,...
            false(size(alp_patterns(:,:,ones(4,1)))),...
            repmat(alp_patterns,[1 1 ndatapoints]));
        alp_patterns_pre = reshape([alp_patterns_pre alp_patterns_pre*0],size(alp_patterns_pre).*[1 1 2]);
        alp_patterns_post = cat(3,...
            repmat(alp_patterns,[1 1 ndatapoints]),...
            false(size(alp_patterns(:,:,ones(4,1)))));
        alp_patterns_post = reshape([alp_patterns_post alp_patterns_post*0],size(alp_patterns_post).*[1 1 2]);
%         switch stim_mode
%             case 'checkerboards'
%                 alp_patterns = cat(3,...
%                     false(size(stim_bg_nospot)),...
%                     alp_patterns,...
%                     check1,...
%                     alp_patterns,...
%                     check2,...
%                     alp_patterns,...
%                     check1,...
%                     alp_patterns,...
%                     check2,...
%                     alp_patterns);
%             case 'bgspots' 
                alp_patterns = repmat(cat(3,...
                    alp_patterns_pre,stim_spot_only,alp_patterns_post,...
                    alp_patterns_pre,stim_bg_nospot,alp_patterns_post,...
                    alp_patterns_pre,stim_bg_nospot+stim_spot_only,alp_patterns_post),[1 1 2]);
%         end
%         a02_load_seq_firefly
%         clear alp_patterns selective_mask stim0_mask stim1_mask stim2_mask stim_bg_nospot stim_spot_only check1 check2 alp_patterns_pre alp_patterns_post
        disp([sequence_mode  ' loaded'])
    case 'structural_hadamard_63c14'
        nstims = 3;
        ndatapoints = 22;
        clc
        blocksize = [63 14];
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,permute(super_mask,[3 2 1]));
        alp_patterns = cat(3,...
            0*alp_patterns(:,:,1),...
            alp_patterns,...
            0*alp_patterns(:,:,ones(2,1)));
        alp_patterns = repmat(alp_patterns,[1 1 nstims]);
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
        a02_load_seq_firefly
        close all
        clear alp_patterns selective_mask stim0_mask stim1_mask stim2_mask
        disp([sequence_mode  ' loaded'])
    case 'original_shorter_63c14'
        nstims = 3;
%         ndatapoints = 22;
        clc
        blocksize = [63 14];
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,permute(super_mask,[3 2 1]));
%         alp_patterns = cat(3,...
%             0*alp_patterns(:,:,1),...
%             alp_patterns,...
%             0*alp_patterns(:,:,ones(2,1)));
%         alp_patterns = repmat(alp_patterns,[1 1 nstims]);
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
%         a02_load_seq_firefly
        close all
%         clear alp_patterns selective_mask stim0_mask stim1_mask stim2_mask
        disp([sequence_mode  ' loaded'])
    case 'shorter_63c14' % the one that gets executed
        blocksize = [63 14];
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,permute(super_mask,[3 2 1]));
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
        alp_patterns = cat(3,alp_patterns*0,alp_patterns);
    case 'used_test_59c8_sim4'
        elementsize = 1;
        blocksizeA = [59 8];
        alp_patternsA = hadamard_patterns_scramble_nopermutation(blocksizeA,elementsize,permute(super_mask,[3 2 1]));
        alp_patternsB = sim_patterns_simpler([],4,1,4);
        alp_patternsB = repmat(alp_patternsB,[1 1 15]);
        alp_patternsC = sim_patterns_simpler([],1,4,4);
        alp_patternsC = repmat(alp_patternsC,[1 1 15]);
        alp_patterns = cat(3,alp_patternsA,reshape([alp_patternsA alp_patternsB alp_patternsA alp_patternsC],size(alp_patternsA).*[1 1 4]));
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
    case 'used_test_59c8_11c3'
        elementsize = 1;
        blocksizeA = [59 8];
        blocksizeB = [11 3];
        alp_patternsA = hadamard_patterns_scramble_nopermutation(blocksizeA,elementsize,permute(super_mask,[3 2 1]));
        alp_patternsB = hadamard_patterns_scramble_nopermutation(blocksizeB,elementsize,permute(super_mask,[3 2 1]));
        alp_patternsB = repmat(alp_patternsB,[1 1 5]);
        alp_patterns = cat(3,alp_patternsA,reshape([alp_patternsA alp_patternsB],size(alp_patternsA).*[1 1 2]));
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
    case 'used_shorter_63c14'
        blocksize = [63 14];
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,permute(super_mask,[3 2 1]));
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
end

%     case 'hadamard_11c3'
%         nstims = 5;
%         ndatapoints = 34;
%         clc
%         blocksize = [11 3];
%         elementsize = 1;
%         alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,permute(super_mask,[3 2 1]));
%         alp_patterns = cat(3,alp_logical_to_btd(permute(super_mask,[2 3 1])),repmat(alp_patterns,[1 1 ndatapoints]));
%         alp_patterns = repmat(alp_patterns,[1 1 nstims]);
%         alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
%         alp_patterns(:,:,4090/nstims*2+1) = alp_logical_to_btd(stim0_mask);
%         alp_patterns(:,:,4090/nstims*3+1) = alp_logical_to_btd(stim1_mask);
%         alp_patterns(:,:,4090/nstims*4+1) = alp_logical_to_btd(stim2_mask);
%         a02_load_seq_firefly
%         close all
%         figure windowstyle docked
%         imshow([selective_mask stim0_mask; stim1_mask stim2_mask],[])
%         title({'selective   stim0','   stim1     stim2'})
%         clear alp_patterns selective_mask stim0_mask stim1_mask stim2_mask
%         disp([sequence_mode  ' loaded'])



end
