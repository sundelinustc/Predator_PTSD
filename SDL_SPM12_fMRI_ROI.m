function SDL_SPM12_fMRI_ROI(SDL)

% Extract ROI data from predetermined ROI (anatomical or functional) into
% Excel datasheet, as well as ANOVA analysis with colorred p value
% Functions are based on Marsbar scripts (http://marsbar.sourceforge.net/faq.html).

% path(genpath(SDL.marsbar_dir),SDL.path); % add path of marsbar

%% Parameters
SDL.roi_path = SDL.ROI_dir;
SDL.fMRI_roi = {
    % Anatomical ROIs, from wfu_pickatlas, HUMAN ATLAS, IBASPM116
%     'L_Amygdala',          'MNI_L_Amygdala_-24_-2_-18_roi.mat';
%     'R_Amygdala',          'MNI_R_Amygdala_27_-1_-19_roi.mat';
%     'L_Hippocampus',       'MNI_L_Hippocampus_-25_-22_-11_roi.mat';
%     'R_Hippocampus',       'MNI_R_Hippocampus_29_-21_-12_roi.mat';
%     'L_ParaHippocampal',   'MNI_L_ParaHippocampal_-21_-17_-22_roi.mat';
%     'R_ParaHippocampal',   'MNI_R_ParaHippocampal_25_-16_-22_roi.mat';
%     
%     'L_Insula',            'MNI_L_Insula_-35_5_2_roi.mat';
%     'R_Insula',            'MNI_R_Insula_39_5_1_roi.mat';
%     'L_IFG',               'MNI_L_IFG_-43_26_5_roi.mat '; % combination of Frontal_Inf_Oper_L, Frontal_Inf_Tri_L and Frontal_Inf_Orib_L
%     'R_IFG',               'MNI_R_IFG_47_26_6_roi.mat'; % combination of Frontal_Inf_Oper_R, Frontal_Inf_Tri_R and Frontal_Inf_Orib_R
%     'L_Rectus',            'MNI_L_Rectus_-5_36_-20_roi.mat';
%     'R_Rectus',            'MNI_R_Rectus_8_34_-19_roi.mat';
%     'L_mOFC',              'MNI_L_mOFC_-5_53_-9_roi.mat'; % Frontal_Mid_Orb_L, the below label (should be Frontal_Med_Orb_L)
%     'R_mOFC',              'MNI_R_mOFC_8_50_-9_roi.mat';  % Frontal_Mid_Orb_R, the below label (should be Frontal_Med_Orb_R)
%     'L_lOFC',              'MNI_L_lOFC_-31_49_-11_roi.mat'; % Frontal_Mid_Orb_L, upper label
%     'R_lOFC',              'MNI_R_lOFC_33_51_-12_roi.mat'; % Frontal_Mid_Orb_R, upper label
%     
%     
%     'L_NAcc',              'MNI_L_NAcc_-13_7_-12_roi.mat'; % from IBASPM71, "nucleus accumbens left"
%     'R_NAcc',              'MNI_R_NAcc_12_9_-11_roi.mat'; % from IBASPM71, "nucleus accumbens right"
%     'L_Caudate',           'MNI_L_Caudate_-12_10_8_roi.mat';
%     'R_Caudate',           'MNI_R_Caudate_14_11_8_roi.mat';
%     'L_Putamen',           'MNI_L_Putamen_-24_3_1_roi.mat';
%     'R_Putamen',           'MNI_R_Putamen_27_4_1_roi.mat';
%     'L_Thalamus',          'MNI_L_Thalamus_-11_-19_7_roi.mat';
%     'R_Thalamus',          'MNI_R_Thalamus_13_-19_7_roi.mat';

%       'L_Precentral',        'MNI_L_Precentral_-39_-7_50_roi.mat';
%       'R_Precentral',        'MNI_R_Precentral_41_-10_51_roi.mat';
%       'L_Postcentral',       'MNI_L_Postcentral_-43_-24_47_roi.mat';
%       'R_Postcentral',       'MNI_R_Postcentral_41_-27_51_roi.mat';
%       'L_SMA',               'MNI_L_SMA_-6_4_60_roi.mat';
%       'R_SMA',               'MNI_R_SMA_8_-1_61_roi.mat';
%       'L_ParaCentral',       'MNI_L_ParaCentral_-8_-27_69_roi.mat';
%       'R_ParaCentral',       'MNI_R_ParaCentral_7_-33_67_roi.mat';
    
    
    % functional ROIs
%     'Raphe',               'MNI_Raphe_-4_-24_-28_roi.mat ';
    % 3 steps of making functional Raphe ROI: (1) fMRI contrast Threat-NonThreat,
    % height-threshold < 0.001, k >= 20 voxels,PTSD < Control,
    % peaking at [-4 -24 -28] mm; (2) a sphere centering at [0 -24 -28] mm
    % (given that raphe is distributed along the midline) with radius = 12
    % mm (3 times of FWHM); (3) r1 & r2.
    
%     'L_Amygdala',          fullfile('Fun','MNI_L_Amygdala_Fun_roi.mat');
%     'R_Amygdala',          fullfile('Fun','MNI_R_Amygdala_Fun_roi.mat');
    % 3 steps of making functional amygdala ROI: 
    % (1) fMRI contrast Threat-NonThreat,
    % height-threshold < 0.025, PTSD < Control,whole brain; 
    % (2) anatomical masks from WFU_PickAtlas for L/R Amygdala; 
    % (3) r1 & r2.
    
    'Andrea_L_IFG',                 fullfile('Andrea','lIFG_handdrawn_-48_30_-3_roi.mat');
    'Andrea_R_IFG',                 fullfile('Andrea','rIFG_handdrawn_r1_r2_roi.mat');
    %     'Andrea_mPFC',                  fullfile('Andrea','mPFC_handdrawn_0_48_4_roi.mat');
    'Andrea_mPFC',                  fullfile('Andrea','mPFC without vmPFC_handdrawn_0_2_48_roi.mat'); % after removing vmPFC
    'Andrea_vmPFC',                 fullfile('Andrea','vmPFC_handdrawn_-0_40_-13_roi.mat');
    %     % These functional masks were from Andrea Gold who is the 1st author of
    %     % the Biological Psychiatry paper "Amygdala-Prefrontal Cortex Functional
    %     % Connectivity During Threat-Induced Anxiety and Goal Distraction"
    };

%       con name        duration (s)
cons = {
    'NonThreat',    30;
    'Threat',       30;
    'Shock',        0.6;
    };

SDL.Cons = {
    'Threat-NonThreat', 'T',0,5,    'con_0004.nii'; % con_0004 is Threat-Nonthreat
    'Threat-None',      'T',0,5,    'con_0002.nii'; % con_0002 is Threat
    'NonThreat-None',      'T',0,5,    'con_0001.nii'; % con_0001 is Nonthreat
    };

%% Percent Signal Change calculation
SDL.psc = SDL.sbjlist(:,1:2); % Subject and Group columns are kept in the new table

for j = 1:size(SDL.fMRI_roi,1) % each ROI
    roi_file = fullfile(SDL.roi_path,SDL.fMRI_roi{j,2}); % e.g. 'C:/Program Files/MATLAB/R2008a/toolbox/addon_for_matlab/marsbar/marsbar-aal-0.2/MNI_Amygdala_L_roi.mat';
    NewColValue = zeros(size(SDL.sbjlist,1),size(cons,1));% the initial value of the % signal change of the contrast of interest
    for m = 1:size(cons,1)
        NewColName{m}  = [SDL.fMRI_roi{j,1},'_',cons{m,1}];
    end
    
    for i = 1:size(SDL.sbjlist,1) % each subject
        spm_name = fullfile(SDL.fMRI_r1st_dir,'Classical',SDL.sbjlist.Subject{i},'SPM.mat');
        
        % Extract % signal change
        D  = mardo(spm_name); % Make marsbar design object
        R  = maroi(roi_file); % Make marsbar ROI object
        Y  = get_marsy(R, D, 'mean'); % Fetch data into marsbar data object
        xCon = get_contrasts(D); % Get contrasts from original design
        E = estimate(D, Y); % Estimate design on ROI data
        E = set_contrasts(E, xCon); % Put contrasts from original design back into design object
        [e_specs, e_names] = event_specs(E); % Get definitions of all events in model
        n_events = size(e_specs, 2);
        
        % Return percent signal esimate for all events in design
        for e_s = 1:n_events
            if strcmp(e_names{e_s},'Shock')
                dur = 0.6;
            else
                dur = 30;
            end
            pct_ev(e_s) = event_signal(E, e_specs(:,e_s), dur);
        end
        
        % Get the mean value of each condition across sessions
         for m = 1:size(cons,1)
            m0 = 0; % the number of the contrast of interest
            for n = 1:size(e_names,2)
                if strcmp(e_names{n},cons{m,1})
                    m0 = m0 + 1;
                    NewColValue(i,m) = NewColValue(i,m) + pct_ev(n);
                end
            end
            if m0 == 0
                NewColValue(i,m) = NaN;
            else
                NewColValue(i,m) = NewColValue(i,m)/m0;
            end
        end
        
        T = array2table(NewColValue,'VariableNames',NewColName); % a new table containing psc of conditions in some ROI across subjects
        [SDL.psc T] % just for display purpose
    end
    SDL.psc = [SDL.psc T];
end

fn = fullfile(SDL.roi_path,'psc.csv'); % excel file can't be used
writetable(SDL.psc,fn);
fprintf('\nPercent Signal Changes saved in: %s\nROI Analyses Completed!!!\n',fn);


%% Percent signal changes between-group comparisons
fprintf('\n=======Begin: Percent signal changes between-group comparisons=======\n');
data = readtable(fullfile(SDL.roi_path,'psc.csv'));

% Beta values, data type, statistical type, statistical test
SDL.test = data.Properties.VariableNames(3:end)'; % the 1st 2 columns are sbj name and group label
% VN = SDL.test(:,1); VN = VN'; % names of variable of interest
% Tmean = varfun(@nanmean,SDL.sbjlist,'InputVariables',VN,'GroupingVariables','Group')
% Tstd  = varfun(@nanstd,SDL.sbjlist, 'InputVariables',VN,'GroupingVariables','Group')


% group name, sbj No. per group
SDL.group = {
    'PTSD',       [26:50]; 
    'Control',    [1:25];
    
    };
% between-group comparisons
% label, g1 and g2 index in SDL.group
SDL.comp = {
    'PTSD vs Control',     1,     2;
    };



fprintf('\n\nTable x. Percent signal changes\nTEST\t%s\t%s\t%s\n',...
    SDL.group{1,1},SDL.group{2,1},SDL.comp{1,1});
fprintf(' \tmean(SD)*\t \t \tt(p)#\n');
for j = 1:size(SDL.test,1)
    % extract data
    fname = SDL.test{j,1}; % e.g. Age
    fdtype = 'num'; % e.g. num = numeric
    jj = find(ismember(data.Properties.VariableNames,fname)); % e.g. 7
    for i = 1:size(data,1)
        fdata = data{i,jj};       
        if strcmp(fdtype,'num')
            if ~isempty(fdata) && isnumeric(fdata)
                SDL.data(i,j) = fdata;
            else
                SDL.data(i,j) = NaN;
            end
        else
        end
    end
    
    % Mean and SD
    for k = 1:size(SDL.group,1)
        % SDL.stat strusture
        % SDL.stat(m1,m2).v1 represents the value out of barackets
        %                .v2                      inside
        % m1 represents SDL.test row
        % m2 represents SDL.group row and more
        % TEST         CONTROL        PTSD        PTSD v CONTROL
        % age          mean(SD)       mean(SD)    t(p)
        % sex          m(f)              m(f)     chi-squ(p)   
 
        fdata = SDL.data(SDL.group{k,2},j);
        fdata = fdata(~isnan(fdata)); % remove NaN
        if     strcmp(fname,'Sex_M1F0')
            SDL.stat(j,k).v1 = sum(fdata==0); % female = 0
            SDL.stat(j,k).v2 = sum(fdata==1); % male   = 1
        else
            SDL.stat(j,k).v1 = mean(fdata);
            SDL.stat(j,k).v2 = std(fdata);
        end
    end
    
    % t and p
    for k = 1:size(SDL.comp,1)
        k0 = size(SDL.group,1); % the t test values are listed following the mean(SD) info
        fg1 = SDL.comp{k,2}; % index of group 1
        fg2 = SDL.comp{k,3}; % index of group 2
        fdata1 = SDL.data(SDL.group{fg1,2},j); % data from group 1
        fdata2 = SDL.data(SDL.group{fg2,2},j); % data from group 2
        fdata1 = fdata1(~isnan(fdata1)); % remove NaN
        fdata2 = fdata2(~isnan(fdata2)); % remove NaN
        if strcmp(fname,'Sex_M1F0') || strcmp(fname,'med_5HT')
            x1 = [repmat('g1',length(fdata1),1); repmat('g2',length(fdata2),1)];
            x2 = [fdata1;fdata2];
            [table,chi2,P] = crosstab(x1,x2);       % chi-square test
            SDL.stat(j,k+k0).v1 = chi2;             % t value
            SDL.stat(j,k+k0).v2 = P;                % p value
        else
            [H,P,CI,STATS] = ttest2(fdata1,fdata2); % 2 sample t test
            SDL.stat(j,k+k0).v1 = STATS.tstat;      % t value
            SDL.stat(j,k+k0).v2 = P;                % p value
        end
    end
    
    % show on the screen
    fprintf('%s\t%1.3f(%1.3f)\t%1.3f(%1.3f)\t%1.3f(%1.3f)\n',...
        fname,...
        SDL.stat(j,1).v1,SDL.stat(j,1).v2,...
        SDL.stat(j,2).v1,SDL.stat(j,2).v2,...
        SDL.stat(j,3).v1,SDL.stat(j,3).v2);
    
end
fprintf('\n=======End: Percent signal changes between-group comparisons=======\n');

end