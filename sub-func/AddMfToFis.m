function out = AddMfToFis(fis, mfParameters, mfTypes,fisVarParameters,fisVarNames,mfNames) 
% this function add the membership function to the fuzzy system

    fuzzyVarParams=fisVarParameters; % get variable parameter
    mf_types_cell= mfTypes;
    mf_parameters=mfParameters;
    fuzzyVarNames = fisVarNames; 
    mf_names = mfNames;
    for i = 1:length(fuzzyVarNames)         % iterate through the number of Variables
        idx = 1;
        fuzzyVar = fuzzyVarParams(i,:); % get variable parameter
        mf_type = mf_types_cell{i};
        sel_var = mf_parameters{i};
        for j = 1:fuzzyVar(2)                    % iterate through the number of MFs 
            if mf_type(j) ==0
                mfparams = sel_var(idx:idx+2);  
                fis = addMF(fis,fuzzyVarNames(i), 'trimf', mfparams, "Name",mf_names{i}(j));
                idx = idx +3; 
            elseif mf_type(j) == 1 
                mfparams = sel_var(idx:idx+3); 
                fis = addMF(fis,fuzzyVarNames(i), 'trapmf', mfparams, "Name",mf_names{i}(j));
                idx = idx +4;
            end
        end
    end
    
    out = fis; 
end