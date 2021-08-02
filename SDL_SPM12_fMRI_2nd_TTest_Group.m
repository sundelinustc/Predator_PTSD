function SDL_SPM12_fMRI_2nd_TTest_Group(SDL)

% SPM12 2nd level model specification, estimation and contrasts setup

spm('Defaults','FMRI');
spm_jobman('initcfg')

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

for j=1:size(Mycon,1)
    
    fprintf('\n2nd-level Betwen-Group Begin: %s\n',Mycon{j,2});
    fdir = fullfile(SDL.fMRI_r2nd_dir,'TTestGroup',Mycon{j,2});
    mkdir(fdir); delete(fullfile(fdir,'SPM.mat'));
    
    % Model specification
    clear matlabbatch;
    matlabbatch{1}.spm.stats.factorial_design.dir = {fdir};
    k1 = 0; k2 = 0;
    for i = 1:size(SDL.sbjlist,1)
        fn = fullfile(SDL.fMRI_r1st_dir,'Classical',SDL.sbjlist.Subject{i},Mycon{j,3});% contrast images
        if strcmp(SDL.sbjlist.Group{i},'PTSD')    % g1 = PTSD
            k1 = k1 + 1;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{k1,1} = fn;
        elseif strcmp(SDL.sbjlist.Group{i},'CONT')% g2 = Control
            k2 = k2 + 1;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{k2,1} = fn;
        end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
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
    fdir = fullfile(SDL.fMRI_r2nd_dir,'TTestGroup',Mycon{j,2});
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(fdir,'SPM.mat')};
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run',matlabbatch);
    
    %% Contrast
    fdir = fullfile(SDL.fMRI_r2nd_dir,'TTestGroup',Mycon{j,2});
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(fdir,'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = [Mycon{j,2},'_g1 ~= g2'];
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = [1 -1];
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = [Mycon{j,2},'_g1 > g2'];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [1 -1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = [Mycon{j,2},'_g1 < g2'];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = [-1 1];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run',matlabbatch);
    fprintf('\n2nd-level Betwen-Group End: %s\n',Mycon{j,2});
end

end