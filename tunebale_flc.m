function out = tunebale_flc(input1,input2,output,mf_types,rules,mf_parameters) 
    %addpath("sub-func","opt-func\");

    %% extract fuzzy variable information 
    
    in1_name = input1.name ;         % input 1 name, fuzzy range and  the number of
    in1_range = input1.range;        % membership functions 
    in1_mfNumber = input1.MfNumber; 
    
    in2_name = input2.name ;         % input 1 name, fuzzy range and  the number of
    in2_range = input2.range;        %  membership functions 
    in2_mfNumber = input2.MfNumber; 
    
    out_name = output.name ;         % output name, fuzzy range and  the number of
    out_range = output.range;        % membership functions 
    out_mfNumber = output.MfNumber; 
    
    fuzzyVarParams = [in1_range, in1_mfNumber;   % Fuzzy variable paramters 
                      in2_range, in2_mfNumber;
                      out_range, out_mfNumber;] ; 
    fuzzyVarNames = [in1_name, in2_name, out_name];  % FIS variable names 

    %% generate a uniformly distributed memebership function
    % generate MFs parameters and MF types
    %mf_parameters = genUniformMfs(fuzzyVarParams,mf_types);
    %mf_parameters  = transformToMfs(fuzzyVarParams, rand_mfs,mf_types);
    %mf_parameters = rand_mfs;
    % generate MFs names 
    mf_names  = genMfNames(mf_types); 
  
    %% Create Fis 
    fis = mamfis('Name','fis');
    fis = addInput(fis,  in1_range, 'Name', in1_name);
    fis = addInput(fis, in2_range, "Name",in2_name); 
    fis = addOutput(fis, out_range,'Name', out_name); 
    
    % add the membership functions to each FIS variable
    fis = AddMfToFis(fis, mf_parameters,mf_types,fuzzyVarParams,fuzzyVarNames,mf_names); 
   
    % add the rule-base of the FIS
    %rule_base = genRuleBase(fuzzyVarParams,rules);
    
    rule_base = rules; 
    fis = addrule(fis, rule_base); 
    out = fis; 
end
% this function generate the names of fuzzy subsets of all variables 
function out = genMfNames(mftypes) 
    mf_names = {["NN"],["NN"],["NN"]};
    fuzzy_vars = ["in1","in2","out"]; 
    cell_size = size(mftypes);
    idx = 1;
    for var = 1: cell_size(2)
        n = size(mftypes{var}); 
                
        for tt =1: n(2)
            mf_names{var}(tt) = fuzzy_vars(var)+"_S"+num2str(tt); 
            idx = idx+1; 
        end 
    end
    out = mf_names; 
end
%this functionadd the membership functions to the FIS 
function out = AddMfToFis(fis, mfParameters, mfTypes,fisVarParameters,fisVarNames,mfNames) 
    
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
        for j = 1:fuzzyVar(3)                    % iterate through the number of MFs 
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