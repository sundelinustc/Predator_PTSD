function SDL_SPM12_fMRI_ROI_corr(SDL)

% (1) Mean and SD (sem) values of gPPI extracted from each ROI
% (2) T-test between 2 groups on the gPPi values
% (3) correlations between gPPI and behavioral performance


% %% Parameters
% % (1) VOI-seed areas
% % (I) VOI name,            (II) VOI image file name
% SDL.VOIList = {
%     'L_Amygdala',          fullfile('Fun','MNI_L_Amy_Fun_roi.nii');
%     'R_Amygdala',          fullfile('Fun','MNI_R_Amy_Fun_roi.nii');
%     % 3 steps of making functional Amygdala ROI: (1) fMRI contrast Threat-NonThreat,
%     % height-threshold < 0.025 (z > 1.96) within L/R Anatomical Amygdala
%     % masks, L Amy peaking at [-20,0,-12] (2 voxels) and [-28,-4,-24] (5
%     % voxels), R Amy peaking at [20,0,-16] (9 voxels); (2) MarsBar r1 | r2 to combine
%     % L Amy 2 cluster into MNI_L_Amy_Fun_roi.mat; (3) MarsBar Export .mat
%     % into .nii images
%     
%     };

% (2) conditions for extracting fMRI
cons = {
    {'Threat'},   {'NonThreat'}, 'T',0,5,    'con_0004.nii'; % con_0004 is Threat-Nonthreat
%     {'Threat'},   {'None'},      'T',0,5,    'con_0002.nii'; % con_0002 is Threat
%     {'NonThreat'},{'None'},      'T',0,5,    'con_0001.nii'; % con_0001 is Nonthreat
    };


% (3) ROI-target areas
SDL.roi = {
    
%     'L_Amygdala',          fullfile('Fun','MNI_L_Amy_Fun_roi.mat');
%     'R_Amygdala',          fullfile('Fun','MNI_R_Amy_Fun_roi.mat');
    % 3 steps of making functional amygdala ROI: 
    % (1) fMRI contrast Threat-NonThreat,
    % height-threshold < 0.025, PTSD < Control,whole brain; 
    % (2) anatomical masks from WFU_PickAtlas for L/R Amygdala; 
    % (3) r1 & r2.
    
%     'Andrea_L_IFG',                 fullfile('Andrea','lIFG_handdrawn_-48_30_-3_roi.mat');
%     'Andrea_R_IFG',                 fullfile('Andrea','rIFG_handdrawn_r1_r2_roi.mat');
%     'Andrea_mPFC',                  fullfile('Andrea','mPFC_handdrawn_0_48_4_roi.mat');
%     'Andrea_mPFC',                  fullfile('Andrea','mPFC without vmPFC_handdrawn_0_2_48_roi.mat'); % after removing vmPFC
    'Andrea_vmPFC',                 fullfile('Andrea','vmPFC_handdrawn_-0_40_-13_roi.mat');    
        % These functional masks were from Andrea Gold who is the 1st author of
        % the Biological Psychiatry paper "Amygdala-Prefrontal Cortex Functional
        % Connectivity During Threat-Induced Anxiety and Goal Distraction"

    
    %     'Ryan_dmPFC',                   fullfile('Ryan','ppi_lAmyg_grpME_dmPFC_-14_58_8_roi.mat');
    %     'Ryan_pgACC',                   fullfile('Ryan','ppi_lAmyg_grpME_pgACC_-1_47_14_roi.mat');
    %     'Ryan_LvmPFC',                  fullfile('Ryan','ppi_lAmyg_grpXage_L_vmPFC_-8_47_-14_roi.mat');
    %     'Ryan_RvmPFC',                  fullfile('Ryan','ppi_lAmyg_grpXage_R_vmPFC_9_44_-16_roi.mat');
    };

% group info
SDL.group = {
    'PTSD';
%     'CONT';
    };

% conditions for correlation analyses
% (I) Behav performance   (II) gPPI values     (III) covariates
SDL.con = {
%     % NonThreat-None
%     'RatePreyCaughtNonThreat',         'NonThreat_minus_None',    {'Age','Sex_M1F0','CTQ','DAST'};
%     'RateAvatarCaughtNonThreat',       'NonThreat_minus_None',    {'Age','Sex_M1F0','CTQ','DAST'};
%     
%     % Threat-None
%     'RatePreyCaughtThreatNonShock',    'Threat_minus_None',       {'Age','Sex_M1F0','CTQ','DAST'};
%     'RateAvatarCaughtThreatNonShock',  'Threat_minus_None',       {'Age','Sex_M1F0','CTQ','DAST'};
    
    % Threat-NonThreat
%     'RatePreyCaughtThreatNonShock_RatePreyCaughtNonThreat',       'Threat_minus_NonThreat',    {'Age','Sex_M1F0','CTQ','DAST'};
    'RateAvatarCaughtThreatNonShock_RateAvatarCaughtNonThreat',   'Threat_minus_NonThreat',    {'Age','Sex_M1F0','CTQ','DAST'};
    };

%% Extract fMRI % signal changes
SDL.psc = SDL.sbjlist(:,{'Subject','Group',...
    'Age','Sex_M1F0','CTQ','DAST',...
    'RatePreyCaughtThreatNonShock_RatePreyCaughtNonThreat','RateAvatarCaughtThreatNonShock_RateAvatarCaughtNonThreat'});
%     'RatePreyCaughtNonThreat','RateAvatarCaughtNonThreat','RatePreyCaughtThreatNonShock','RateAvatarCaughtThreatNonShock'});


% for j = 1:size(SDL.VOIList,1)             % for each VOI
%     fVOI = SDL.VOIList{j,1};
    
    for k = 1:size(cons,1)                % for each contrast of interest
        fCon = [char(cons{k,1}),'_minus_',char(cons{k,2})];
        
        for m = 1:size(SDL.roi,1)         % for each ROI
            fROI = SDL.roi{m,1};
            NewColName  = [fCon,'_',fROI]; % column name in the data table
            NewColValue = zeros(size(SDL.sbjlist,1),1);% the initial value of the gPPI values of the contrast of interest
            
            for i = 1:size(SDL.sbjlist,1) % for each subject
                roi_files = fullfile(SDL.ROI_dir,SDL.roi{m,2});
%                 fdPPI = fullfile(SDL.psc_r1st_dir,'Classical',SDL.sbjlist.Subject{i},['PPI_',fVOI]); % directory containing the contrast image
%                 fnPPI = ['con_PPI_',fCon,'_',fVOI,'.nii']; % name of contrast image
                P = fullfile(SDL.fMRI_r1st_dir,'Classical',SDL.sbjlist.Subject{i},cons{k,6});% contrast images
                rois = maroi('load_cell', roi_files);  % make maroi ROI objects
                mY = get_marsy(rois{:}, P, 'mean');  % extract data into marsy data object
                y = summary_data(mY); % get summary time course(s)
                fprintf('%d/%d subjects completed\n',i,size(SDL.sbjlist,1));
                
                NewColValue(i,1) = y;
                T = array2table(NewColValue,'VariableNames',{NewColName}); % a new table containing psc of conditions in some ROI across subjects
                %                 [SDL.psc T]; % for display purpose only
            end
            SDL.psc = [SDL.psc T];
            
        end
    end
% end
SDL.psc % for display purpose
fn = fullfile(SDL.fMRI_r2nd_dir,'fMRI_BetaValue_corr.csv'); % excel file can't be used
writetable(SDL.psc,fn);
fprintf('\ngPPI Beta Values saved in: %s\n',fn);


%% T tests: mean(SD) and t(p) values
fprintf('\n\nTable 1. fMRI Beta Values\ngPPIs\tPTSD [mean(SD)]\tCONT [mean(SD)]\t t-value(p-value)\n');
for iCon = 1:size(cons,1) % for each of conditions to be compared, i.e. Threat-none or NonThreat-none
%     for iVOI = 1:size(SDL.VOIList,1) % for each seed area
        
        plist = [];
        for iROI = 1:size(SDL.roi,1) % for each target area
            
            % name of the output
            fname = [cons{iCon,1}{1},'_minus_',cons{iCon,2}{1},'_',SDL.roi{iROI,1}];
            iPTSD = find(ismember(SDL.psc.Group,'PTSD')); % find the index of PTSD group members
            iCONT = find(ismember(SDL.psc.Group,'CONT')); % find the index of CONT group members
            
            mean_PTSD = mean(SDL.psc{iPTSD,fname});std_PTSD = std(SDL.psc{iPTSD,fname});
            mean_CONT = mean(SDL.psc{iCONT,fname});std_CONT = std(SDL.psc{iCONT,fname});
            [H,P,CI,STATS] = ttest2(SDL.psc{iPTSD,fname},SDL.psc{iCONT,fname}); % 2 sample t test
            plist = [plist P];
            fprintf('%s\t%1.3f(%1.3f)\t%1.3f(%1.3f)\t%1.3f(%1.3f)\n',fname,...
                mean_PTSD,std_PTSD,mean_CONT,std_CONT,...
                STATS.tstat,P);
        end
        
        % FDR correction of p values --only for Threat_minus_NonThreat
        [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(plist,.05,'pdep','yes');     % seed=L or R Amy
        fprintf('\nFDR correction\ttadj_p: '); adj_p
        fprintf('\n');
%     end
end


%% correlation analyses: R(p) and between-group t(p) values
fprintf('\n\nTable 2. Correlation between fMRI beta values and behavioral performance.\nBehav\tgPPIs\tPTSD [R(p)]\tCONT [R(p)]\t between-group t-value(p-value)\n');
for iCon = 1:size(SDL.con,1) % for each of conditions to do correlation analyses
%     for iVOI = 1:size(SDL.VOIList,1) % for each seed area
        plist = [];
        for iROI = 1:size(SDL.roi,1) % for each target area
            
            % name of the output
            fname1 = SDL.con{iCon,1};
            fname2 = [SDL.con{iCon,2},'_',SDL.roi{iROI,1}];
            fname3 = SDL.con{iCon,3};
            iPTSD = find(ismember(SDL.psc.Group,'PTSD')); % find the index of PTSD group members
            iCONT = find(ismember(SDL.psc.Group,'CONT')); % find the index of CONT group members
            
            [R,p] = partialcorr(SDL.psc{iPTSD,{fname1,fname2,fname3{:}}}); R_PTSD=R(1,2); p_PTSD=p(1,2); % correlation between col 1&2 is at matrix(1,2)
            [R,p] = partialcorr(SDL.psc{iCONT,{fname1,fname2,fname3{:}}}); R_CONT=R(1,2); p_CONT=p(1,2); % correlation between col 1&2 is at matrix(1,2)
            
            [z,p] = SDL_CorrCmp2(R_PTSD,R_CONT,length(iPTSD),length(iCONT));
            plist = [plist p];
            fprintf('%s\t%s\t%1.3f(%1.3f)\t%1.3f(%1.3f)\t%1.3f(%1.3f)\n',fname1,fname2,...
                R_PTSD,p_PTSD,R_CONT,p_CONT,...
                z,p);
        end
        
        % FDR correction of p values --only for Threat_minus_NonThreat
        [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(plist,.05,'pdep','yes');     % seed=L or R Amy
        fprintf('\nFDR correction\ttadj_p: '); adj_p
        fprintf('\n');
%     end
end


fprintf('\n\nAll analyses completed!!!\n\n');

end