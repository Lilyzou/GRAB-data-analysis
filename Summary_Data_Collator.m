%% Notes
% This file is designed to scrape data from the behavior file and raw
% fluorescence from the load file, putting it into structures. 
% Please load the load file manually before running
% this script. Run the script for each session you want to scrape data
% from.


%%
n = 1;
continue_read = true;
trials = b{1,1}.trials;
while (continue_read)
    R(n).f = tmp; %#ok<*SAGROW> % Get raw F
    R(n).trialType = zeroes(4,length(trials));
    for i=1:length(trials) %%%% not sure what to do here 
        switch trials{1,i}.trialType
            case 1
                R(n).trialType(1,i) = 1;
            case 2
                R(n).trialType(2,i) = 1;
            case 3
                R(n).trialType(3,i) = 1;
            case 4
                R(n).trialType(4,i) = 1;
        end
    end
    
    
    
    toContinue = input('Continue? Y/N');
    if (toContinue == 'N' || toContinue == 'n')
        continue_read = false;
    else
        n = n + 1;
    end
        
end