function SDL_SPM12_fMRI_2nd_Flex_Session_Group(SDL)

% SPM12 2nd level model specification, estimation and contrasts setup

spm('Defaults','FMRI');
spm_jobman('initcfg')

Mycon = {
    'NonThreat',        {'con_0001.nii','con_0004.nii','con_0007.nii','con_0010.nii','con_0013.nii'};
    'Threat',           {'con_0002.nii','con_0005.nii','con_0008.nii','con_0011.nii','con_0014.nii'};
    'Shock',            {'con_0003.nii','con_0006.nii','con_0009.nii','con_0012.nii','con_0015.nii'};
    'Threat-NonThreat', {'con_0016.nii','con_0017.nii','con_0018.nii','con_0019.nii','con_0020.nii'},...
    };

for j=1:size(Mycon,1)
    
    fprintf('\n2nd-level Betwen-Group Begin: %s\n',Mycon{j,1});
    fdir = fullfile(SDL.fMRI_r2nd_dir,'SessionInter','BetweenGroup',['Flex_',Mycon{j,1}]);
    mkdir(fdir); delete(fullfile(fdir,'SPM.mat'));
    
    % Model specification
    clear matlabbatch;
    matlabbatch{1}.spm.stats.factorial_design.dir = {fdir};
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'subject';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
    nc = 5;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Session'; % within-subject factor, 5 levels: 1-5 sessions
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
    ng = 2;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'Group'; % between-group factor, 2 levels: PTSD vs Control
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
    MEc = [1:nc]-mean(1:nc); % (main effect of condition, here: [-2 -1 0 1 2]) 
    MEg = [1 -1]; % (main effect of group: Group 1 > Group 2)
    
    n1 = 0; n2 = 0; % subjects' number in each group
    for i = 1:size(SDL.sbjlist,1)
        sdir = fullfile(SDL.fMRI_r1st_dir,'SessionInter',SDL.sbjlist.Subject{i});
        if strcmp(Mycon{j,1},'NonThreat')
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).scans = {
                fullfile(sdir,'con_0001.nii,1');
                fullfile(sdir,'con_0004.nii,1');
                fullfile(sdir,'con_0007.nii,1');
                fullfile(sdir,'con_0010.nii,1');
                fullfile(sdir,'con_0013.nii,1');
                };
        elseif strcmp(Mycon{j,1},'Threat')
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).scans = {
                fullfile(sdir,'con_0002.nii,1');
                fullfile(sdir,'con_0005.nii,1');
                fullfile(sdir,'con_0008.nii,1');
                fullfile(sdir,'con_0011.nii,1');
                fullfile(sdir,'con_0014.nii,1');
                };
        elseif strcmp(Mycon{j,1},'Shock')
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).scans = {
                fullfile(sdir,'con_0003.nii,1');
                fullfile(sdir,'con_0006.nii,1');
                fullfile(sdir,'con_0009.nii,1');
                fullfile(sdir,'con_0012.nii,1');
                fullfile(sdir,'con_0015.nii,1');
                };
        elseif strcmp(Mycon{j,1},'Threat-NonThreat')
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).scans = {
                fullfile(sdir,'con_0016.nii,1');
                fullfile(sdir,'con_0017.nii,1');
                fullfile(sdir,'con_0018.nii,1');
                fullfile(sdir,'con_0019.nii,1');
                fullfile(sdir,'con_0020.nii,1');
                };
        else
        end
        
        if strcmp(SDL.sbjlist.Group{i},'PTSD')
            n1 = n1 + 1;
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).conds = [
                1 1
                1 2
                1 3
                1 4
                1 5
                ];
        elseif strcmp(SDL.sbjlist.Group{i},'CONT')
            n2 = n2 + 1;
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(i).conds = [
                2 1
                2 2
                2 3
                2 4
                2 5
                ];
        else
        end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.inter.fnums = [2
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
    fdir = fullfile(SDL.fMRI_r2nd_dir,'SessionInter','BetweenGroup',['Flex_',Mycon{j,1}]);
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(fdir,'SPM.mat')};
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run',matlabbatch);
    
    %% Contrast
    MyCon = {
%         'T', 'PTSD-CONT,sess01',                [1 0 0 0 0 -1 0 0 0 0 ones(1,n1)/n1 -1*ones(1,n2)/n2]; 
%         'T', 'CONT-PTSD,sess01',                [-1 0 0 0 0 1 0 0 0 0 -1*ones(1,n1)/n1 ones(1,n2)/n2]; 
%         
%         'T', 'PTSD-CONT,sess02',                [0 1 0 0 0 0 -1 0 0 0 ones(1,n1)/n1 -1*ones(1,n2)/n2];
%         'T', 'CONT-PTSD,sess02',                [0 -1 0 0 0 0 1 0 0 0 -1*ones(1,n1)/n1 ones(1,n2)/n2];
%         
%         'T', 'PTSD-CONT,sess03',                [0 0 1 0 0 0 0 -1 0 0 ones(1,n1)/n1 -1*ones(1,n2)/n2];
%         'T', 'CONT-PTSD,sess03',                [0 0 -1 0 0 0 0 1 0 0 -1*ones(1,n1)/n1 ones(1,n2)/n2];
%         
%         'T', 'PTSD-CONT,sess04',                [0 0 0 1 0 0 0 0 -1 0 ones(1,n1)/n1 -1*ones(1,n2)/n2];
%         'T', 'CONT-PTSD,sess04',                [0 0 0 -1 0 0 0 0 1 0 -1*ones(1,n1)/n1 ones(1,n2)/n2];
%         
%         'T', 'PTSD-CONT,sess05',                [0 0 0 0 1 0 0 0 0 -1 ones(1,n1)/n1 -1*ones(1,n2)/n2];
%         'T', 'CONT-PTSD,sess05',                [0 0 0 0 -1 0 0 0 0 1 -1*ones(1,n1)/n1 ones(1,n2)/n2];
% 
%         'T', 'PTSD-CONT,Early',                 [1 1 0 0 0 -1 -1 0 0 0 2*ones(1,n1)/n1 -2*ones(1,n2)/n2];
%         'T', 'CONT-PTSD,Early',                 [-1 -1 0 0 0 1 1 0 0 0 -2*ones(1,n1)/n1 2*ones(1,n2)/n2];
%         
%         'T', 'PTSD-CONT,Late',                  [0 0 0 1 1 0 0 0 -1 -1 2*ones(1,n1)/n1 -2*ones(1,n2)/n2];
%         'T', 'CONT-PTSD,Late',                  [0 0 0 -1 -1 0 0 0 1 1 -2*ones(1,n1)/n1 2*ones(1,n2)/n2];
        
        'T', 'PTSD>CONT,Late-Early',            [-1 -1 0 1 1 1 1 0 -1 -1 zeros(1,n1) zeros(1,n2)];
        'T', 'CONT>PTSD,Late-Early',            [1 1 0 -1 -1 -1 -1 0 1 1 zeros(1,n1) zeros(1,n2)];
        };
    fdir = fullfile(SDL.fMRI_r2nd_dir,'SessionInter','BetweenGroup',['Flex_',Mycon{j,1}]);
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(fdir,'SPM.mat')};
    for ic = 1:size(MyCon,1) 
        if     strcmp(MyCon{ic,1}, 'F')
            matlabbatch{1}.spm.stats.con.consess{ic}.fcon.name    = MyCon{ic,2};
            matlabbatch{1}.spm.stats.con.consess{ic}.fcon.convec  = MyCon{ic,3}; % e.g.{[ones(1,n1)/n1 -ones(1,n2)/n2 MEg zeros(1,nc) ones(1,nc)/nc -ones(1,nc)/nc]};
            matlabbatch{1}.spm.stats.con.consess{ic}.fcon.sessrep = 'none';
        elseif strcmp(MyCon{ic,1}, 'T')
            matlabbatch{1}.spm.stats.con.consess{ic}.tcon.name    = MyCon{ic,2};
            matlabbatch{1}.spm.stats.con.consess{ic}.tcon.convec  = MyCon{ic,3};
            matlabbatch{1}.spm.stats.con.consess{ic}.tcon.sessrep = 'none';
        else
            fprintf('\nWrong contrast type!!!\n');
        end
    end
    
    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run',matlabbatch);
    fprintf('\n2nd-level Betwen-Group End: %s\n',Mycon{j,1});
end

end