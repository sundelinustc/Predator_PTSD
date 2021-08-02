function SDL_SPM12_fMRI_2nd_Flex_Session_EachGroup(SDL)

% SPM12 2nd level model specification, estimation and contrasts setup

spm('Defaults','FMRI');
spm_jobman('initcfg')

Group = {
    'CONT'
    'PTSD'
    };

Mycon0 = {
    'T',  'NonThreat_sess01',          [zeros(1,9*0) [1 0 0 0 0 0 0 0 0] zeros(1,9*4)]; % 1
    'T',  'Threat_sess01',             [zeros(1,9*0) [0 1 0 0 0 0 0 0 0] zeros(1,9*4)]; % 2
    'T',  'Shock_sess01',              [zeros(1,9*0) [0 0 1 0 0 0 0 0 0] zeros(1,9*4)]; % 3
    
    'T',  'NonThreat_sess02',          [zeros(1,9*1) [1 0 0 0 0 0 0 0 0] zeros(1,9*3)]; % 4
    'T',  'Threat_sess02',             [zeros(1,9*1) [0 1 0 0 0 0 0 0 0] zeros(1,9*3)]; % 5
    'T',  'Shock_sess02',              [zeros(1,9*1) [0 0 1 0 0 0 0 0 0] zeros(1,9*3)]; % 6
    
    'T',  'NonThreat_sess03',          [zeros(1,9*2) [1 0 0 0 0 0 0 0 0] zeros(1,9*2)]; % 7
    'T',  'Threat_sess03',             [zeros(1,9*2) [0 1 0 0 0 0 0 0 0] zeros(1,9*2)]; % 8
    'T',  'Shock_sess03',              [zeros(1,9*2) [0 0 1 0 0 0 0 0 0] zeros(1,9*2)]; % 9
    
    'T',  'NonThreat_sess04',          [zeros(1,9*3) [1 0 0 0 0 0 0 0 0] zeros(1,9*1)]; % 10
    'T',  'Threat_sess04',             [zeros(1,9*3) [0 1 0 0 0 0 0 0 0] zeros(1,9*1)]; % 11
    'T',  'Shock_sess04',              [zeros(1,9*3) [0 0 1 0 0 0 0 0 0] zeros(1,9*1)]; % 12
    
    'T',  'NonThreat_sess05',          [zeros(1,9*4) [1 0 0 0 0 0 0 0 0] zeros(1,9*0)]; % 13
    'T',  'Threat_sess05',             [zeros(1,9*4) [0 1 0 0 0 0 0 0 0] zeros(1,9*0)]; % 14
    'T',  'Shock_sess05',              [zeros(1,9*4) [0 0 1 0 0 0 0 0 0] zeros(1,9*0)]; % 15
    };

% Cov = {
%     'Age';
%     'Sex_M1F0';
%     'CTQ';
%     'DAST';
%     };

tlist = SDL.sbjlist;
for ip = 1:size(Group,1)
    SDL.sbjlist = tlist;
    ind = strncmp(SDL.sbjlist.Group,Group{ip,1},4); % look for the rows matching the group name
    SDL.sbjlist(~ind,:)=[]; % delete the rows non-matching the group name
    
    
    fprintf('\n2nd-level Begin: Group=%s\n',Group{ip,1});
    fdir = fullfile(SDL.fMRI_r2nd_dir,'SessionInter',['Flex_',Group{ip,1}]);
    mkdir(fdir); delete(fullfile(fdir,'SPM.mat'));
    
    % Model specification
    clear matlabbatch;
    matlabbatch{1}.spm.stats.factorial_design.dir = {fdir};
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'subject';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Stress'; % 2 levels: Threat vs NonThreat
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'Session'; % 5 levels: 1-5 sessions
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
    for i = 1:size(SDL.sbjlist,1)
        sdir = fullfile(SDL.fMRI_r1st_dir,'SessionInter',SDL.sbjlist.Subject{i});
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).scans = {
            fullfile(sdir,'con_0001.nii,1');
            fullfile(sdir,'con_0002.nii,1');
            fullfile(sdir,'con_0004.nii,1');
            fullfile(sdir,'con_0005.nii,1');
            fullfile(sdir,'con_0007.nii,1');
            fullfile(sdir,'con_0008.nii,1');
            fullfile(sdir,'con_0010.nii,1');
            fullfile(sdir,'con_0011.nii,1');
            fullfile(sdir,'con_0013.nii,1');
            fullfile(sdir,'con_0014.nii,1');
            };
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).conds = [
            1 1 % NonThreat,sess01
            1 2 % Threat,sess01
            2 1 % NonThreat,sess02
            2 2 % Threat,sess02
            3 1 % NonThreat,sess03
            3 2 % Threat,sess03
            4 1 % NonThreat,sess04
            4 2 % Threat,sess04
            5 1 % NonThreat,sess05
            5 2 % Threat,sess05
            ];
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [2
        3];
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    spm_jobman('run',matlabbatch);
    
    
    %% Estimation
    fdir = fullfile(SDL.fMRI_r2nd_dir,'SessionInter',['Flex_',Group{ip,1}]);
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(fdir,'SPM.mat')};
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run',matlabbatch);
    
    %% Contrast
    MyCon = {
        'T', '(Late>Early),NonThreat',                [-1 0 -1 0 0 0 1 0 1 0]; %  SPMT_0001.nii
        'T', '(Late<Early),NonThreat',                [1 0 1 0 0 0 -1 0 -1 0]; %  SPMT_0002.nii
        'T', '(Late>Early),Threat',                   [0 -1 0 -1 0 0 0 1 0 1]; %  SPMT_0003.nii
        'T', '(Late<Early),Threat',                   [0 1 0 1 0 0 0 -1 0 -1]; %  SPMT_0004.nii
        'T', '(Late>Early)*(Threat-NonThreat)',       [1 -1 1 -1 0 0 -1 1 -1 1];% SPMT_0005.nii
        'T', '(Late<Early)*(Threat-NonThreat)',       [-1 1 -1 1 0 0 1 -1 1 -1];% SPMT_0006.nii
        }; 
    
    fdir = fullfile(SDL.fMRI_r2nd_dir,'SessionInter',['Flex_',Group{ip,1}]);
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(fdir,'SPM.mat')};
    for j = 1:size(MyCon,1) 
        if     strcmp(MyCon{j,1}, 'F')
            matlabbatch{1}.spm.stats.con.consess{j}.fcon.name    = MyCon{j,2};
            matlabbatch{1}.spm.stats.con.consess{j}.fcon.convec  = MyCon{j,3}; % e.g.{[ones(1,n1)/n1 -ones(1,n2)/n2 MEg zeros(1,nc) ones(1,nc)/nc -ones(1,nc)/nc]};
            matlabbatch{1}.spm.stats.con.consess{j}.fcon.sessrep = 'none';
        elseif strcmp(MyCon{j,1}, 'T')
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.name    = MyCon{j,2};
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec  = MyCon{j,3};
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
        else
            fprintf('\nWrong contrast type!!!\n');
        end
    end
    
    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run',matlabbatch);
    
end

end