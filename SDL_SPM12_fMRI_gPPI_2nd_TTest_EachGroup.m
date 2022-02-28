function SDL_SPM12_fMRI_gPPI_2nd_TTest_EachGroup(SDL)

% The script does 2nd level analysis for gPPI model
% Notice! Other toolbox such as eeglab or marsbar may disrupt the SPM working, so, DO remove them when using SPM

pack; % increase memory
spm('Defaults','FMRI');
spm_jobman('initcfg');

SDL.Group = {
    'CONT'
    'PTSD'
    };

% (1)VOI name, (2)mask images
SDL.VOIList = {    
  % Anatomical ROIs, from wfu_pickatlas, HUMAN ATLAS, IBASPM116
%     'L_Amygdala',          'MNI_L_Amygdala.nii';
%     'R_Amygdala',          'MNI_R_Amygdala.nii';

%     'L_Caudate',          'MNI_L_Caudate.nii';
%     'R_Caudate',          'MNI_R_Caudate.nii';

%       'L_ACC',          'MNI_L_ACC.nii';
%       'R_ACC',          'MNI_R_ACC.nii';

% 'L_Hippocampus',          'MNI_L_Hippocampus.nii';
% 'R_Hippocampus',          'MNI_R_Hippocampus.nii';
% 
% 'L_ParaHippocampal',      'MNI_L_ParaHippocampal.nii';
% 'R_ParaHippocampal',      'MNI_R_ParaHippocampal.nii';

  
%     % functional ROIs
%     'Raphe',               'MNI_Raphe.nii';
%     % 3 steps of making functional Raphe ROI: (1) fMRI contrast Threat-NonThreat,
%     % height-threshold < 0.001, k >= 20 voxels,PTSD < Control,
%     % peaking at [-4 -24 -28] mm; (2) a sphere centering at [0 -24 -28] mm
%     % (given that raphe is distributed along the midline) with radius = 12
%     % mm (3 times of FWHM); (3) MarsBar r1 & r2.

    'L_Amygdala',          fullfile('Fun','MNI_L_Amy_Fun_roi.nii');
    'R_Amygdala',          fullfile('Fun','MNI_R_Amy_Fun_roi.nii');
    % 3 steps of making functional Amygdala ROI: (1) fMRI contrast Threat-NonThreat,
    % height-threshold < 0.025 (z > 1.96) within L/R Anatomical Amygdala
    % masks, L Amy peaking at [-20,0,-12] (2 voxels) and [-28,-4,-24] (5
    % voxels), R Amy peaking at [20,0,-16] (9 voxels); (2) MarsBar r1 | r2 to combine
    % L Amy 2 cluster into MNI_L_Amy_Fun_roi.mat; (3) MarsBar Export .mat
    % into .nii images
    
    };

SDL.Cons = {
    {'Threat'},   {'NonThreat'}, 'T',0,5;
    {'Threat'},   {'None'},      'T',0,5;
    {'NonThreat'},{'None'},      'T',0,5;
    };

Cov = {
    'Age';
    'Sex_M1F0';
    'CTQ';
    'DAST';
    };

tlist = SDL.sbjlist;
for ip = 1:size(SDL.Group,1)
    ind = strncmp(SDL.sbjlist.Group,SDL.Group{ip,1},4); % look for the rows matching the group name
    SDL.sbjlist = tlist;
    SDL.sbjlist(~ind,:)=[]; % delete the rows non-matching the group name
    
    
    for j = 1:size(SDL.VOIList,1)
        for k = 1:size(SDL.Cons,1)
            if     size(SDL.Cons{k,1},2) == 1 % only 1 condition, for simple main effect
                conname = [SDL.Cons{k,1}{1},'_minus_',SDL.Cons{k,2}{1},'_']; % always left minus right
            elseif size(SDL.Cons{k,1},2) == 2 % 2 conditions, for interaction effect
                conname = [SDL.Cons{k,1}{1},'_',SDL.Cons{k,1}{2},'_minus_',SDL.Cons{k,2}{1},'_',SDL.Cons{k,2}{2},'_']; % always left minus right
            else
            end
            %% Model specification
            fdir = fullfile(SDL.gPPI_r2nd_dir,['TTest_',SDL.Group{ip,1}],conname,SDL.VOIList{j,1});
            mkdir(fdir); delete(fullfile(fdir,'SPM.mat'));
            
            clear matlabbatch
            matlabbatch{1}.spm.stats.factorial_design.dir = {fdir};
            k1 = 0;
            for i = 1:size(SDL.sbjlist,1) % for each subject
                fn = fullfile(SDL.gPPI_r1st_dir,'Classical',SDL.sbjlist.Subject{i},...
                    ['PPI_',SDL.VOIList{j,1}],['con_PPI_',conname,SDL.VOIList{j,1},'.nii']);% contrast images
                if exist(fn,'file') == 2 % check if a file exists -- some gPPI files do not exist, e.g.20140905_18841, L Amygdala
                    k1 = k1 + 1;
                    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{k1,1} = fn;
                else
                    fprintf('Warning: No such a file: %s\n',fn);
                end
            end
            for k1 = 1:size(Cov,1)
                fname = Cov{k1,1};
                jj = find(ismember(SDL.sbjlist.Properties.VariableNames,fname)); % e.g. 7
                matlabbatch{1}.spm.stats.factorial_design.cov(k1).c = SDL.sbjlist{:,jj};
                matlabbatch{1}.spm.stats.factorial_design.cov(k1).cname = Cov{k1};
                matlabbatch{1}.spm.stats.factorial_design.cov(k1).iCFI = 1;
                matlabbatch{1}.spm.stats.factorial_design.cov(k1).iCC = 1;
            end
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            spm_jobman('run',matlabbatch);
            
            %% Model Estimation
            clear matlabbatch;
            matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(fdir,'SPM.mat')};
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            spm_jobman('run',matlabbatch);
            
            %% Contrast
            clear matlabbatch;
            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(fdir,'SPM.mat')};
            % matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = {''};
            % matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = {1};
            % matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = ['PPI_',conname,'> 0'];                   % spmT_0001.img
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1];
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = ['PPI_',conname,'< 0'];                   % spmT_0002.img
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1];
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            spm_jobman('run',matlabbatch);
            
        end
    end
    
end

end