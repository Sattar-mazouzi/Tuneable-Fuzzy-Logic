function [state, options, optchanged] = captureBestCosts(options, state, flag)
    persistent bestCosts meanCosts
    
    if strcmp(flag, 'init')
        bestCosts = [];
        meanCosts = [];
    end
    
    if strcmp(flag, 'iter')
        bestCosts(end + 1) = min(state.Score);
        meanCosts(end + 1) = mean(state.Score);
    end
    
    if strcmp(flag, 'done')
        assignin('base', 'bestCostsAcrossGenerations', bestCosts);
        assignin('base', 'meanCostsAcrossGenerations', meanCosts);
    end
    
    optchanged = false;
end