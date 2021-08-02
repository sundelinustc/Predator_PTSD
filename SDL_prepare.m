function SDL_prepare(SDL)

% prepare files for preprocessing
% -- copyfiles from carr to my desktop, to accelerate process
% -- unzip .nii.gz files, for SPM

%% parameters
%   raw filename,           corresponding target folder name
flist = {
    'run01.nii.gz',         'FunImg';
    'run02.nii.gz',         'S2_FunImg';
    'run03.nii.gz',         'S3_FunImg';
    'run04.nii.gz',         'S4_FunImg';
    'run05.nii.gz',         'S5_FunImg';
    'T1_brain.nii.gz',      'T1Img';
    };

% %% copy, unzip and delete files
% for i=1:size(SDL.sbjlist,1) % for each subject
%     for j = 1:size(flist,1) % for each file
%         fn1 = fullfile(SDL.fMRI_raw_dir, SDL.sbjlist.Subject{i},'rawdata',flist{j,1});
%         fn2 = fullfile(SDL.fMRI_orig_dir,SDL.sbjlist.Subject{i});mkdir(fn2); % mkdir first
%         fn2 = fullfile(SDL.fMRI_orig_dir,SDL.sbjlist.Subject{i},          flist{j,1});
%         % copyfile
%         copyfile(fn1,fn2);
%         fprintf('Copyfile:%s -->\n\t%s\n',fn1,fn2);
%         % gunzip
%         gunzip(fn2);
%         fprintf('Unzip:%s\n',fn2);
%         % delete .nii.gz file
%         delete(fn2);
%         fprintf('Delete:%s\n\n',fn2);
%     end
% end

% %% 4D to 3D saved in Preprocess directory
% for i=1:size(SDL.sbjlist,1) % for each subject
%     for j = 1:size(flist,1)-1 % for each functional file
%         clear matlabbatch;
%         fn1 = fullfile(SDL.fMRI_orig_dir,SDL.sbjlist.Subject{i},flist{j,1}(1:end-3));
%         fn2 = fullfile(SDL.fMRI_prep_dir,flist{j,2},SDL.sbjlist.Subject{i});mkdir(fn2);
%         matlabbatch{1}.spm.util.split.vol = {[fn1,',1']};
%         matlabbatch{1}.spm.util.split.outdir = {fn2};
%         spm_jobman('run',matlabbatch);
%         fprintf('4D->3D:%s -->\n\t%s\n',fn1,fn2);
%     end
% end

% %% copy T1 images to Preprocess directory
% for i=1:size(SDL.sbjlist,1) % for each subject
%     for j = size(flist,1) % for the T1 file -- the last file in fList
%         fn1 = fullfile(SDL.fMRI_orig_dir,SDL.sbjlist.Subject{i},flist{j,1}(1:end-3));
%         fn2 = fullfile(SDL.fMRI_prep_dir,flist{j,2},SDL.sbjlist.Subject{i});mkdir(fn2);
%         copyfile(fn1,fn2);
%         fprintf('Copyfile:%s -->\n\t%s\n',fn1,fn2);
%     end
% end


% %% L-R flip
% % given that Courtney has converted both functional and T1 images into LAP
% % orientation for FSL, while SPM uses RAP orientation
% for i=1:size(SDL.sbjlist,1) % for each subject
%     for j = 1:size(flist,1) % for each file (including both functional and T1 images)
%         % L-R flip
%         clear matlabbatch;
%         fn1 = fullfile(SDL.fMRI_prep_dir,flist{j,2},SDL.sbjlist.Subject{i});
%         fn2 = cellstr(spm_select('FPList',fn1,'.*.nii'));
%         matlabbatch{1}.spm.util.reorient.srcfiles = fn2;
%         matlabbatch{1}.spm.util.reorient.transform.transM = [
%             -1 0 0 0
%             0 1 0 0
%             0 0 1 0
%             0 0 0 1
%             ]; % matrix for L-R flip
%         matlabbatch{1}.spm.util.reorient.prefix = 'f'; % prefix for flipping
%         spm_jobman('run',matlabbatch);
%         fprintf('L-R flip:%s\n',fn1);
%         
%         % Movefile to new folder
%         fn3 = fullfile(SDL.fMRI_prep_dir,[flist{j,2},'f'],SDL.sbjlist.Subject{i});
%         mkdir(fn3);
%         fn2 = cellstr(spm_select('FPList',fn1,'^f.*.nii'));
%         for k = 1:size(fn2,1)
%             [fp,fn,fe] = fileparts(fn2{k,1});
%             movefile(fn2{k,1},fullfile(fn3,[fn,fe]));
%             fprintf('Movefile:%s ->\n\t%s\n',fn2{k,1},fullfile(fn3,[fn,fe]));
%         end
%     end
% end


%% Copy behavioral files
fn = fullfile(SDL.fMRI_prep_dir,'Behav'); % directory containing behavioral data
mkdir(fn); fprintf('Creat Directory: %s\n',fn);
for i=1:size(SDL.sbjlist,1) % for each subject
    fn1 = fullfile('/Volumes/data/Predator.01/Analysis/fsl2/',SDL.sbjlist.Subject{i},'behavioral');
    fn2 = fullfile(fn,SDL.sbjlist.Subject{i});
    copyfile(fn1,fn2);
    fprintf('Copyfile: %s\n----> %s\n',fn1,fn2);  
end

end