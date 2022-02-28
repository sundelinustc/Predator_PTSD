function SDL_SPM12_fMRI_gPPI_1st(SDL)

% The script sets gPPI model and does 1st level analysis in choice stage for all participant
% Notice! Other toolbox such as eeglab or marsbar may disrupt the SPM working, so, DO remove them when using SPM

% pack; % increase memory
spm('Defaults','FMRI');
spm_jobman('initcfg');
% addpath(fullfile(SDL.Batch_dir,'misc')); % the function 'SDL_create_sphere_image.m' is needed

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

for i = 1:2 %size(SDL.sbjlist,1)
    workdir = fullfile(SDL.gPPI_r1st_dir,'Classical',SDL.sbjlist.Subject{i});
    cd(workdir);
    
    for j = 1:size(SDL.VOIList,1)
        fprintf('\n1st-level gPPI Begin: %dth, sbj=%s, VOI=%s\n',i, SDL.sbjlist.Subject{i},SDL.VOIList{j,1});
        
        %% Generate PPPI structure
        try
            rmdir(fullfile(workdir,['PPI_' SDL.VOIList{j,1}]),'s');% delete the exusted VOI folder
            fprintf('\nFolder Removed: %s\n',fullfile(workdir,['PPI_' SDL.VOIList{j,1}]));
        catch
            fprintf('\nFolder NON-EXIST: %s\n',fullfile(workdir,['PPI_' SDL.VOIList{j,1}]));
        end
        
        fn1 = fullfile(SDL.ROI_dir,[SDL.VOIList{j,2}]);
        fn2 = fullfile(SDL.fMRI_r2nd_dir,'TTest','Shock','mask.nii,1');
        
  
        % logic AND the VOI image and whole brain mask image
        % to avoid the VOI image going beyond the extent for all subjects
        % This is only for R STG
        clear matlabbatch;
        matlabbatch{1}.spm.util.imcalc.input = {
            fn1; % e.g. 'H:\SDL\Delin_SH_2014\Analysis\fMRI\1st\choice\sbj01\VOI_R Str_21_-3_0_mask.nii,1'
            fn2; % e.g. 'H:\SDL\Delin_SH_2014\Analysis\fMRI\2nd\ANOVA\choice\mask.img,1'
            };
        matlabbatch{1}.spm.util.imcalc.output = [SDL.VOIList{j,1},'_mask_andWholeBrainMask']; % e.g. 'VOI_R Str_21_-3_0_mask_andWholeBrainMask';
        matlabbatch{1}.spm.util.imcalc.outdir = {workdir};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1 & i2'; % logic AND
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch);
        clear matlabbatch;
        fprintf('VOI and whole brain mask made: %s\n',[SDL.VOIList{j,1},'_mask_andWHoleBrainMask.nii']);
        
        
        % set the PPPI parameter structure
        P = [];
        P.subject = SDL.VOIList{j,1};
        P.directory = workdir;
        P.VOI = fullfile(workdir,[SDL.VOIList{j,1},'_mask_andWholeBrainMask.nii']);
        P.Region = SDL.VOIList{j,1};
        P.Estimate = 1;
        P.contrast = 0;
        P.extract = 'eig'; % eig=eigenvariate, mean=mean values
        P.Tasks = {'1'  'NonThreat'  'Threat'  'Shock'};
        P.Weights = [];
        P.analysis = 'psy';
        P.method = 'cond'; % cond=gPPI method, trad=traditional SPM method
        P.CompContrasts = 1;
        P.Weighted = 0;
        
        for k = 1:size(SDL.Cons,1)
            P.Contrasts(k).left      = SDL.Cons{k,1};
            P.Contrasts(k).right     = SDL.Cons{k,2};
            P.Contrasts(k).STAT      = SDL.Cons{k,3};
            P.Contrasts(k).Weighted  = SDL.Cons{k,4};
            P.Contrasts(k).MinEvents = SDL.Cons{k,5};
            if     size(SDL.Cons{k,1},2) == 1 % only 1 condition, for simple main effect
                P.Contrasts(k).name = [SDL.Cons{k,1}{1} '_minus_' SDL.Cons{k,2}{1}]; % always left minus right
            elseif size(SDL.Cons{k,1},2) == 2 % 2 conditions, for interaction effect
                P.Contrasts(k).name = [SDL.Cons{k,1}{1},'_',SDL.Cons{k,1}{2}     '_minus_'     SDL.Cons{k,2}{1},'_',SDL.Cons{k,2}{2}]; % always left minus right
            else
            end
        end
        
        % run PPPI
        P
        PPPI(P);
        clear P;
        fprintf('\n1st-level gPPI End: %dth, sbj=%s, VOI=%s\n',i, SDL.sbjlist.Subject{i},SDL.VOIList{j,1});
    end
end


end