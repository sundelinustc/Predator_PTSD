function SDL_SPM12_fMRI_1st(SDL)

% SPM12 1st level model specification, estimation and contrasts setup

% the folders/files containing images (1st col), behav data (2nd col) and realign parameter (3rd col) for each run
SDL.sess = {
    'FunImgARWS',      'run1',      'rp_afrun01_00001.txt';
    'S2_FunImgARWS',   'run2',   'S2_rp_afrun02_00001.txt';
    'S3_FunImgARWS',   'run3',   'S3_rp_afrun03_00001.txt';
    'S4_FunImgARWS',   'run4',   'S4_rp_afrun04_00001.txt';
    'S5_FunImgARWS',   'run5',   'S5_rp_afrun05_00001.txt';
    };

% condition name,               regressor files included
SDL.con = {
    'NonThreat'                 ,{'NonThreat'};
    'Threat'                    ,{'ThreatWithoutShock','ThreatWithShock'};
    'Shock'                     ,{'Shock'};
    };

dir1 = {
    'Classical'; % the folder containg the 1st-level analyses
%     'TEST';
    }; % for test

for i=1:size(SDL.sbjlist,1) % for each subject
    fprintf('\n1st-level Begin: %s\n',SDL.sbjlist.Subject{i});
    
    %% 1) Model Specification
    fn = fullfile(SDL.fMRI_r1st_dir,dir1{1},SDL.sbjlist.Subject{i});
    mkdir(fn); fprintf('Mkdir: %s\n',fn);
    delete(fullfile(fn,'SPM.mat'));
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_spec.dir = {fn};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = SDL.TR; % 2
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = SDL.nslices; % 34, Num of time bins per slice
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = SDL.nslices/2; %, middle of time bins per slice
    
    for j = 1:size(SDL.sess,1) % for each session
        fImg  = fullfile(SDL.fMRI_prep_dir,dir1{1},SDL.sess{j,1},SDL.sbjlist.Subject{i}); % folders containing images
        fBehv = fullfile(SDL.fMRI_prep_dir,dir1{1},'Behav',SDL.sbjlist.Subject{i},SDL.sess{j,2}); % folders containing behavioral data
        matlabbatch{1}.spm.stats.fmri_spec.sess(j).scans = cellstr(spm_select('FPList',fImg,'^sw.*\.nii'));
        
        for k = 1:size(SDL.con,1) % for each condition
            fdata = [];
            for m = 1:size(SDL.con{k,2},2)
                fdata = [fdata;load(fullfile(fBehv,[SDL.con{k,2}{m},'.stf']))];
            end
            fdata = sortrows(fdata,1); % sort data acording to onset time
            if strcmp(SDL.con{k,1},'NoThreat') || strcmp(SDL.con{k,1},'NonThreat')
                fdata(:,1) = fdata(:,1) - 2; % onset time should be the cue onset time 2s ahead (acording to Andrea's suggestion)
                fdata(:,2) = fdata(:,2) + 2; % block duration for Threat and NonThreat Blocks are 32s (according to Andrea's suggestion)
            end
            fhead = fullfile(SDL.fMRI_prep_dir,dir1{1},'RealignParameter',SDL.sbjlist.Subject{i},SDL.sess{j,3}); % headmotion relignment
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).cond(k).name = SDL.con{k,1}; % e.g. 'NonThreat'
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).cond(k).onset = fdata(:,1);
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).cond(k).duration = fdata(1,2); % one value is OK
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).cond(k).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).cond(k).orth = 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi_reg = {fhead};
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).hpf = 128;
        end
        
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    spm_jobman('run',matlabbatch);
    
    
    %% 2) Estimation
    fn = fullfile(SDL.fMRI_r1st_dir,dir1{1},SDL.sbjlist.Subject{i});
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(fn,'SPM.mat')};
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run',matlabbatch);
    
    %% 3) Contrast
    % Basic condition
    Mycon0 = {
        'NonThreat',                 repmat([1 0 0 0 0 0 0 0 0],1,5);
        'Threat',                    repmat([0 1 0 0 0 0 0 0 0],1,5);
        'Shock',                     repmat([0 0 1 0 0 0 0 0 0],1,5);
        };
    % Contrasts of Interest
    Mycon = {
        'T',  'NonThreat',                         Mycon0{1,2};              %1
        'T',  'Threat',                            Mycon0{2,2};              %2
        'T',  'Shock',                             Mycon0{3,2};              %3
        'T',  'Threat-NonThreat',                  Mycon0{2,2}-Mycon0{1,2};  %4
        'T',  'Shock-NonThreat',                   Mycon0{3,2}-Mycon0{1,2};  %5
        'T',  'Shock-Threat',                      Mycon0{3,2}-Mycon0{2,2};  %6
        };
    fn = fullfile(SDL.fMRI_r1st_dir,dir1{1},SDL.sbjlist.Subject{i});
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(fn,'SPM.mat')};
    for j = 1:size(Mycon,1)
        if strcmp(Mycon{j,1},'F')
            matlabbatch{1}.spm.stats.con.consess{j}.fcon.name = Mycon{j,2};
            matlabbatch{1}.spm.stats.con.consess{j}.fcon.convec = Mycon{j,3};
            matlabbatch{1}.spm.stats.con.consess{j}.fcon.sessrep = 'none';
        else
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.name = Mycon{j,2};
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = Mycon{j,3};
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
        end
    end
    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run',matlabbatch);
    
    fprintf('\n1st-level Completed: %s\n',SDL.sbjlist.Subject{i});
end



end