% Function to check that the user input matches available options
%
% INPUTS
% phrase   = the phrase that will be displayed by input (string)
% options  = the valid options for user responses (cell array of strings)
%
% OUTPUTS  
% s        = the string input by the user
%
% CHANGE LOG
% Created 2020/05/26 RML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = inputChecker(phrase,options)

    u = input(phrase);
    
    goodInput = sum(strcmp(u,options));
    
    if goodInput > 0
        s = u;
    else
        optionString = '';
        for kk = 1:length(options)
            optionString = [optionString ', ''' options{kk} ''''];
        end
        disp(' ')
        disp(['Please choose from the following options for your input:' optionString])
        s = inputChecker(phrase,options);
    end

end