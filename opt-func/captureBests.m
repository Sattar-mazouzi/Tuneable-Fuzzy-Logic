function [state, options, optchanged] = captureBests(options, state, flag)
    persistent bestCosts meanCosts bestPopulations
    
    if strcmp(flag, 'init')
        bestCosts = [];
        meanCosts = [];
        bestPopulations = [];
    end
    
    if strcmp(flag, 'iter')
        bestCosts(end + 1) = min(state.Score);
        meanCosts(end + 1) = mean(state.Score);
        % Save the best population for the current iteration
        [~, bestIdx] = min(state.Score); % Find the index of the best score
        bestPopulations(end + 1, :) = state.Population(bestIdx, :); % Save the best individual
    end
    
    if strcmp(flag, 'done')
        assignin('base', 'bestCostsAcrossGenerations', bestCosts);
        assignin('base', 'meanCostsAcrossGenerations', meanCosts);
        assignin('base', 'bestPopulationsAcrossGenerations', bestPopulations);
    end
    
    optchanged = false;
end