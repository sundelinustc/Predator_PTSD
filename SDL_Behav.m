function SDL_Behav(SDL)


%% Demographic variables between-group comparisons
fprintf('\n=======Begin: Demographic information=======\n');

% demographic variable, data type, statistical type, statistical test
SDL.test = {
    % Demographic information
    'Age',                'num',   'mean(SD)',       't test';
    'Sex_M1F0',           'num',   'sum(0)/sum(1)',  'chi square';
    'CAPS_C',             'num',   'mean(SD)',       't test';
    'CAPS_L',             'num',   'mean(SD)',       't test';
    'BDI',                'num',   'mean(SD)',       't test';
    'CTQ',                'num',   'mean(SD)',       't test';
    'DAST',               'num',   'mean(SD)',       't test';
    'STAI_state',         'num',   'mean(SD)',       't test';
    'STAI_trait',         'num',   'mean(SD)',       't test';
    'DRRI_Sec1',          'num',   'mean(SD)',       't test';
    'DRRI_Sec2',          'num',   'mean(SD)',       't test';
    'DRRI_Sec3',          'num',   'mean(SD)',       't test';
    'DRRI_Sec4',          'num',   'mean(SD)',       't test';

    % Behavioral performance
    'Difficulty',         'num',   'mean(SD)',       't test';
%     'Score_Run1',         'num',   'mean(SD)',       't test';
%     'Score_Run2',         'num',   'mean(SD)',       't test';
%     'Score_Run3',         'num',   'mean(SD)',       't test';
%     'Score_Run4',         'num',   'mean(SD)',       't test';
%     'Score_Run5',         'num',   'mean(SD)',       't test';
%     'SbjCaught_Run1',     'num',   'mean(SD)',       't test';
%     'SbjCaught_Run2',     'num',   'mean(SD)',       't test';
%     'SbjCaught_Run3',     'num',   'mean(SD)',       't test';
%     'SbjCaught_Run4',     'num',   'mean(SD)',       't test';
%     'SbjCaught_Run5',     'num',   'mean(SD)',       't test';
%     'PreyCaught_Run1',    'num',   'mean(SD)',       't test';
%     'PreyCaught_Run2',    'num',   'mean(SD)',       't test';
%     'PreyCaught_Run3',    'num',   'mean(SD)',       't test';
%     'PreyCaught_Run4',    'num',   'mean(SD)',       't test';
%     'PreyCaught_Run5',    'num',   'mean(SD)',       't test';
    'RatePreyCaughtNonThreat',   'num',   'mean(SD)',       't test';
    'RateAvatarCaughtNonThreat', 'num',   'mean(SD)',       't test';
    
    'RatePreyCaughtThreatNonShock','num',   'mean(SD)',       't test';
    'RateAvatarCaughtThreatNonShock','num',   'mean(SD)',       't test';
    
    'RatePreyCaughtThreat',      'num',   'mean(SD)',       't test'; % after excluding the time points during shock and 12s post the shock
    'RateAvatarCaughtThreat',    'num',   'mean(SD)',       't test'; % after excluding the time points during shock and 12s post the shock
   
    };
% VN = SDL.test(:,1); VN = VN'; % names of variable of interest
% Tmean = varfun(@nanmean,SDL.sbjlist,'InputVariables',VN,'GroupingVariables','Group')
% Tstd  = varfun(@nanstd,SDL.sbjlist, 'InputVariables',VN,'GroupingVariables','Group')


% group name, sbj No. per group
SDL.group = {
    'PTSD',     [26:50];
    'Control',  [1:25];
    };
% between-group comparisons
% label, g1 and g2 index in SDL.group
SDL.comp = {
    'PTSD vs Control',     1,     2;
    };



fprintf('\n\nTable 1. Demographic Information\nTEST\t%s\t%s\t%s\n',...
    SDL.group{1,1},SDL.group{2,1},SDL.comp{1,1});
fprintf(' \tmean(SD)*\t \t \tt(p)#\n');
for j = 1:size(SDL.test,1)
    % extract data
    fname = SDL.test{j,1}; % e.g. Age
    fdtype = SDL.test{j,2}; % e.g. num = numeric
    jj = find(ismember(SDL.sbjlist.Properties.VariableNames,fname)); % e.g. 7
    for i = 1:size(SDL.sbjlist,1)
        fdata = SDL.sbjlist{i,jj};       
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
fprintf('Note: *, values out/in brackets are number of females/males for "sex". # statistical values are from Chi-Square tests for "sex"\n');
fprintf('\n=======End: Demographic information=======\n');





end