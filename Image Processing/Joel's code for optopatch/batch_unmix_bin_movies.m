clear
baseDir = 'C:\Joel\Data\';
expDir = {'2014-07-22 3d Round 1 Drugs (c)\Take 10,RB9dCorr,000uM riluzole',...
    '2014-07-22 3d Round 1 Drugs (c)\Take 8-RB9d,010uM riluzole',...
    '2014-07-22 3d Round 1 Drugs (c)\Take 9-RB9d,100uM riluzole',...
%     '2014-07-22 3d Round 1 Drugs (c)\Take 4-39bCorr,000uM riluzole',...
%     '2014-07-22 3d Round 1 Drugs (c)\Take 5-39bCorr,010uM riluzole',...
%     '2014-07-22 3d Round 1 Drugs (c)\Take 6-39bCorr,100uM riluzole',...
    }

% '2014-03-13 Round 1 ALS (2d diff)\Take 1-2.5,12dpp,9dpi,Tyr',...
%     '2014-03-13 Round 1 ALS (2d diff)\Take 2-39B,12dpp,9dpi,Tyr',...
%     '2014-03-13 Round 1 ALS (2d diff)\Take 3-19f,12dpp,9dpi,Tyr',...
%     '2014-03-13 Round 1 ALS (2d diff)\Take 4-Hus,12dpp,9dpi,Tyr'
% '2014-02-14 Rat primary shRNA (rapa)\Take 1-Con,16dpp,14dpi2V,shCon,tyr+bl',...
%     '2014-02-14 Rat primary shRNA (rapa)\Take 2-shTSC,15dpp,14dpi2V,shTSC,tyr+bl',...
%     '2014-02-14 Rat primary shRNA (rapa)\Take 3-shTSC+rapa_16dpp,14dpi,2V,shTSC,tyr+bl'};
% '2014-02-11 Rat primary shRNA\Take 1-Con,12dpp,11dpi2VshRNA,Tyr+blockers',...
%     '2014-02-11 Rat primary shRNA\Take 2-shTSC,12dpp,11dpi2VshRNA,Tyr+blockers'}
% '2014-02-05 TSC2 hES (no glia)\Take 1-WT,54dpd,18dpi2V,Tyr',...
%     '2014-02-05 TSC2 hES (no glia)\Take 2-HET,54dpd,18dpi2V,Tyr'};
% '2014-02-05 Rat primrary (shRNA)\Take 1-Con,16dpp,14dpiCon,12dpi2V,Tyr+bloc',...
%     '2014-02-05 Rat primrary (shRNA)\Take 2-KD,16dpp,14dpiTSC,12dpi2V,Tyr+bloc',...
% '2014-02-03 shRNA rat primary\Take 1-shTSC,12dpp,10dpSH,8dp2V,Tyr+blockers',...
%     '2014-02-03 shRNA rat primary\Take 2-shCON,12dpp,10dpSH,8dp2V,Tyr+blockers',...
%     '2014-02-03 TSC2 hES (glia)\Take 1-WT,52dpp,21dpi2V,Tyr+block',...
%     '2014-02-03 TSC2 hES (glia)\Take 2-HET,52dpp,21dpi2V,Tyr+block'};
%      '2014-01-24 retigabine again\Take 1-17dpp,13dpi,2V-Blockers'
%     '2014-01-26 TSC2 hES\Take 1-WT,44dpd,16dpi,2V',...
%     '2014-01-26 TSC2 hES\Take 3-HET,44dpd,16dpi,2V',...
%     '2014-01-27 shRNA rat primary\Take 1-con,18dpp,17dpiCON,13dpi2V,Tyr+bl',...
%     '2014-01-27 shRNA rat primary\Take 2-sh,18dpp,17dpiSH,13dpi2V,Tyr+bl',...
%     '2014-01-27 TSC2 hES (glia)\Take 1-WT,45dpd,17dpi2V,Tyr+blockers',...
%     '2014-01-27 TSC2 hES (glia)\Take 2-Het,45dpd,17dpi2V,Tyr+blockers',...
%     '2014-01-27 TSC2 hES (glia)\Take 3-KO,45dpd,17dpi2V,Tyr+blockers'};

% '2014-01-23 shRNA cells\Take 1-Con_16dpp,14dpshRNA,12dp2V',...
%     '2014-01-23 shRNA cells\Take 2-TSC_16dpp,14dpshRNA,12dp2V,+rapa'};%,...

nexpDir = length(expDir);
for iD = 1:nexpDir
    runDir = [baseDir expDir{iD}];
    dname = dir(runDir);
    ndir = length(dname)-2;
    
    matchdir = '\d\d';
    
    for ii = 1:ndir
        testDir = ~isempty(regexp(dname(ii+2).name(1:2),matchdir));
        if dname(ii+2).isdir == 1 && testDir
            ['working on file ' num2str(ii) ' of ' num2str(ndir)...
                ' in directory ' num2str(iD) ' of ' num2str(nexpDir)]
            ncells = 6;
            nr = 104;
            nc = 118;
            cycleLen = 525;
            currDir = pwd;
            
            t1 = tic;
            [avgImgRed, GFP, traces, cellimgs, filters, rawTrace] = ...
                unmix_bin_optopatch_movies(runDir,ii,ncells,nr,nc,cycleLen);
            disp(['entire file read took ' num2str(toc(t1)) ' s']);
            
            cd(currDir);
        end
    end
end