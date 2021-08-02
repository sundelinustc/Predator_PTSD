function SDL_SPM12_fMRI_2nd_TTest_EachGroup(SDL)

% SPM12 2nd level model specification, estimation and contrasts setup

spm('Defaults','FMRI');
spm_jobman('initcfg')

Group = {
    'CONT'
    'PTSD'
    };

Mycon = {
    'T',  'NonThreat',                           'con_0001.nii';
    'T',  'Threat',                              'con_0002.nii';
    'T',  'Shock',                               'con_0003.nii';
    'T',  'Threat-NonThreat',                    'con_0004.nii';
    'T',  'Shock-NonThreat',                     'con_0005.nii';
    'T',  'Shock-Threat',                        'con_0006.nii';
    };

Cov = {
    'Age';
    'Sex_M1F0';
    'CTQ';
    'DAST';
    };

tlist = SDL.sbjlist;
for ip = 1:size(Group,1)
    SDL.sbjlist = tlist;
    ind = strncmp(SDL.sbjlist.Group,Group{ip,1},4); % look for the rows matching the group name
    SDL.sbjlist(~ind,:)=[]; % delete the rows non-matching the group name
    
    for j=1:size(Mycon,1)
        
        fprintf('\n2nd-level Begin: Group=%s,con=%s\n',Group{ip,1},Mycon{j,2});
        fdir = fullfile(SDL.fMRI_r2nd_dir,['TTest_',Group{ip,1}],Mycon{j,2});
        mkdir(fdir); delete(fullfile(fdir,'SPM.mat'));
        
        % Model specification
        clear matlabbatch;
        matlabbatch{1}.spm.stats.factorial_design.dir = {fdir};
        for i = 1:size(SDL.sbjlist,1)
            fn = fullfile(SDL.fMRI_r1st_dir,'Classical',SDL.sbjlist.Subject{i},Mycon{j,3});% contrast images
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{i,1} = fn;
        end
        for k = 1:size(Cov,1)
            fname = Cov{k,1};
            jj = find(ismember(SDL.sbjlist.Properties.VariableNames,fname)); % e.g. 7
            matlabbatch{1}.spm.stats.factorial_design.cov(k).c = SDL.sbjlist{:,jj};
            matlabbatch{1}.spm.stats.factorial_design.cov(k).cname = Cov{k};
            matlabbatch{1}.spm.stats.factorial_design.cov(k).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(k).iCC = 1;
        end
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        spm_jobman('run',matlabbatch);
        
        %% Estimation
        fdir = fullfile(SDL.fMRI_r2nd_dir,['TTest_',Group{ip,1}],Mycon{j,2});
        clear matlabbatch;
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(fdir,'SPM.mat')};
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('run',matlabbatch);
        
        %% Contrast
        fdir = fullfile(SDL.fMRI_r2nd_dir,['TTest_',Group{ip,1}],Mycon{j,2});
        clear matlabbatch;
        matlabbatch{1}.spm.stats.con.spmmat = {fullfile(fdir,'SPM.mat')};
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = Mycon{j,2};
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = 1;
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = [Mycon{j,2},' > 0'];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = 1;
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = [Mycon{j,2},' < 0'];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = -1;
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.delete = 1;
        spm_jobman('run',matlabbatch);
        
    end
end

end