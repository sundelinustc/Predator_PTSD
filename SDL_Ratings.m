% Load and organize ratings from each participant.
% There are ratings for 8 questons at the end of the task.
%   1	For the *BLUE* mazes, did you focus more on avoiding predator (=1) or catching prey (=7)?
%   2	How anxious were you when you entered *BLUE* mazes? 1=not at all to 9=highly anxious
%   3	During the *BLUE* mazes, how anxious were you? 1=not at all to 9=highly anxious
%   4	During the *BLUE* mazes, how much did you dread being chased by predator? 1=not at all to 9=high dread
%   5	For the *RED* mazes, did you focus more on avoiding predator (=1) or catching prey (=7)?
%   6	How anxious were you when you entered *RED* mazes? 1=not at all to 9=highly anxious
%   7	During the *RED* mazes, how anxious were you? 1=not at all to 9=highly anxious
%   8	During *RED* mazes, how much did you dread being chased by predator? 1=not at all to 9=high dread
% The output should be a Table, each row represents a subject, each column a question answer

T = SDL.sbjlist(:,{'Subject','Group'}); % table to containing the outputs

% loop to access the ratings
for i = 1:size(T,1) % per subject
    if i==39 || i==40 || i==41
    else
        fn  = fullfile(SDL.P_dir,'Data','GPF',T.Subject{i},['GPF_',T.Subject{i}(end-4:end),'_10.out']); % filename of the ratings
        fid = fopen(fn,'r');
        A   = textscan(fid,'%s'); % loading all text info
        if i == 47 % sbj 20008 only has 5 ratings
            for j = 1:5 % per rating (totally 5)
                T{i,j+2} = str2num(A{1,1}{j+1}(end)); % line 2 to 9 are ratings
            end
        else
            for j = 1:8 % per rating (totally 8)
                T{i,j+2} = str2num(A{1,1}{j+1}(end)); % line 2 to 9 are ratings
            end
        end
        fclose(fid);
    end  
end
T.Properties.VariableNames = {'Subject', 'Group', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8'}; % change the column names