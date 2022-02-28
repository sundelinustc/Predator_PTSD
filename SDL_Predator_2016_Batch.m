% function SDL_Predator_2016_Batch

% For fMRI study of Predator, Duke University,2016-07-12
% by Dr. Delin Sun

%% Parameters
% pack;
% clear SDL;

% [nhpath, nhfile, nhext] = fileparts(which('spm'));
% SDL.SPM_dir =     nhpath; % e.g. '/home/neuropsy2/Downloads/MatlabToolBox/spm12/',  SPM12, version 6225, server: Neuropsy2
% [nhpath, nhfile, nhext] = fileparts(which('marsbar'));
% SDL.marsbar_dir = nhpath;
% SDL.path = path;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project paths
SDL.P_dir =          fullfile('/Volumes/data/Predator.01'); % carr 
% SDL.P1_dir =         fullfile('/Users/ds366/Desktop/temp/Predator'); % my own desktop
SDL.P1_dir =         fullfile('/Volumes/Morey/Lab/Delin/Projects/Predator/SPM'); % munin
SDL.Batch_dir =      fullfile(SDL.P1_dir,'Scripts');
SDL.Raw_dir =        fullfile(SDL.P_dir, 'Analysis','fsl2');
SDL.Original_dir =   fullfile(SDL.P1_dir,'Original');
SDL.Preprocess_dir = fullfile(SDL.P1_dir,'Preprocess');
SDL.Analysis_dir =   fullfile(SDL.P1_dir,'Analysis');
SDL.ROI_dir =        fullfile(SDL.P1_dir,'ROI');

%% Sbjstcs' information
SDL.sbjlist = readtable(fullfile(SDL.Batch_dir,'predator_subjects.xlsx'),'Sheet','Sheet1');
% ind = find(strcmp(SDL.sbjlist.Include,'N')); % look for excluded subjects
ind = find(SDL.sbjlist.Include_Y1N0==0);
if ~isempty(ind)
    SDL.sbjlist(ind,:) = []; % remove the excluded subjects
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% fMRI Task part %%%%%%%%%%%
%% fMRI parameters
SDL.fMRI_raw_dir  = SDL.Raw_dir; % .nii.gz files are stored in carr
SDL.fMRI_orig_dir = SDL.Original_dir;          % NIFTI 
SDL.fMRI_prep_dir = fullfile(SDL.Preprocess_dir,'fMRI');
SDL.fMRI_r1st_dir = fullfile(SDL.Analysis_dir ,'fMRI','1st');
SDL.fMRI_r2nd_dir = fullfile(SDL.Analysis_dir ,'fMRI','2nd');
% SDL.fMRI_roi = cellstr(spm_select('FPList',SDL.marsbar_aal_dir,'.*roi.mat')); % all MNI aal anatomical ROI

SDL.TR = 2; % TR = 2s
SDL.nslices = 34; % slices number = 34
SDL.TA = SDL.TR-SDL.TR/SDL.nslices; % TA = TR - TR/nslices
SDL.slice_order = [1:2:33,2:2:34]; % slices acquisation order, for slice timing
SDL.refslice = 33; % reference slice, i.e. the middle one of the SDL.slice_order

cd(fullfile(SDL.Batch_dir,'fMRI'));

% Preprocess
% SDL_prepare(SDL); % copy and unzip files

% Behavioral data analyses
% SDL_Behav(SDL); 
SDL_Ratings(SDL);

% 1st level analyses
% SDL_SPM12_fMRI_1st(SDL);      % Classical analyses
% SDL_SPM12_fMRI_1st_Session(SDL); % Investigating Threat/NonThreat/Shock x Sesion

% 2nd level analyses
% SDL_SPM12_fMRI_2nd_TTest(SDL); % all subjects togather
% SDL_SPM12_fMRI_2nd_TTest_EachGroup(SDL); % test within each group
% SDL_SPM12_fMRI_2nd_TTest_Group(SDL); % Between-group differences
% SDL_SPM12_fMRI_2nd_Corr_EachGroup(SDL); % correlation between fMRI signals and performance

% SDL_SPM12_fMRI_2nd_Flex_Session_EachGroup(SDL); % test within each group, (Threat-NonThreat) x Sesion
% SDL_SPM12_fMRI_2nd_Flex_Session_Group(SDL); % Between-group differences, (Threat-NonThreat) x Sesion

% ROI analyses
% SDL_SPM12_fMRI_ROI(SDL); % using % signal change, should be replaced by beta values in SDL_SPM12_fMRI_ROI_corr(SDL)
% SDL_SPM12_fMRI_ROI_Session(SDL); % Threat/NonThreat/Shock x Sesion
% SDL_SPM12_fMRI_ROI_corr(SDL); % beta values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% gPPI part %%%%%%%%%%%
% gPPI parameters
% gPPI paths
SDL.gPPI_r1st_dir = fullfile(SDL.Analysis_dir ,'gPPI','1st');
SDL.gPPI_r2nd_dir = fullfile(SDL.Analysis_dir ,'gPPI','2nd');
% gPPI analysis
% cd(fullfile(SDL.Batch_dir,'gPPI'));
% SDL_SPM12_fMRI_1st(SDL); % re-run in a new directory to differ from fMRI gPPI results
% SDL_SPM12_fMRI_gPPI_1st(SDL);
% SDL_SPM12_fMRI_gPPI_2nd_TTest_EachGroup(SDL); % test within each group
% SDL_SPM12_fMRI_gPPI_2nd_TTest_Group(SDL); % Between-group differences
% SDL_SPM12_fMRI_gPPI_2nd_Corr_EachGroup(SDL); % correlation between gPPI (seed:L/R Amygdala) and performance
% SDL_SPM12_gPPI_ROI(SDL); % extracting PPI beta coefficient from contrast images
% SDL_SPM12_gPPI_ROI_corr(SDL); % correlations between gPPI and behavioral performance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Resting State part %%%%%%%%%%%
% RESTING parameters
% REST paths
SDL.nslices = 30; % slices number = 40
SDL.TA = 2-2/30; % TA = TR - TR\nslices
SDL.slice_order = [1:2:29,2:2:30]; % slices acquisation order in SH, for slice timing
SDL.refslice = 29; % reference slice, i.e. the middle one of the SDL.slice_order
SDL.REST_prep_dir = fullfile(SDL.Preprocess_dir,'REST'); % e.g. 'H:\SDL\Delin_SH_2014\Preprocess\REST\';
SDL.REST_ana_dir =  fullfile(SDL.Analysis_dir,'REST');

% REST Analysis
% cd(fullfile(SDL.Batch_dir,'REST'));
% SDL_SPM8_fMRI_Rest_prepare(SDL);
% SDL_SPM8_fMRI_Rest_Task_Choice(SDL);
% SDL_SPM12_REST_prepare1214(SDL); % prepare data for SH2012 & 2014 analyses


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% cat12 part %%%%%%
% SBM parameters
% SBM paths
SDL.SBM_prep_dir =  fullfile(SDL.Preprocess_dir,'SBM');
SDL.SBM_ana_dir =   fullfile(SDL.Analysis_dir,  'SBM');
SDL.VBM_ana_dir =   fullfile(SDL.Analysis_dir,  'VBM');

% SBM Analysis
% cd(fullfile(SDL.Batch_dir,'SBM'));
% SDL_SPM12_preprocessing_cat12(SDL);
% SDL_SPM12_SBM_2nd_TTest_Group(SDL); % SBM between-group comparisons
% SDL_SPM12_VBM_2nd_TTest_Group(SDL); % VBM between-group comparisons
% SDL_SPM12_VBM_2nd_Corr_EachGroup(SDL); % correlation analyses
% SDL_SPM12_VBM_ROI_corr(SDL);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% DTI part %%%%%%
% cd([SDL.Batch_dir    'DTI']);
% SDL_SPM8_fMRI_DTI_prepare(SDL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Mediation part %%%%%%
%% Mediation parameters
SDL.Mediation_dir =         fullfile(SDL.Analysis_dir,  'Mediation');
SDL.RobustRegression_dir =  fullfile(SDL.Analysis_dir,  'RobustRegression');

%% Mediation Analysis
% cd(fullfile(SDL.Batch_dir,'Mediation'));


% SDL_RobustRegression(SDL);
% SDL_Mediation_Analysis(SDL);

%% End
% cd(SDL.Batch_dir);