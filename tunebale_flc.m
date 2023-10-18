function out = tunebale_flc(input1,input2,output,mf_types,rules,rand_mfs ) 
    addpath("sub-func","opt-func\");

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
    [mf_parameters, mf_types_cell] = genUniformMfs(fuzzyVarParams,mf_types); 
    mf_parameters  = transformToMfs(fuzzyVarParams, rand_mfs,mf_types_cell);
    % generate MFs names 
    mf_names  = genMfNames(mf_types_cell); 
  
    %% Create Fis 
    fis = mamfis('Name','fis');
    fis = addInput(fis,  in1_range, 'Name', in1_name);
    fis = addInput(fis, in2_range, "Name",in2_name); 
    fis = addOutput(fis, out_range,'Name', out_name); 
    
    % add the membership functions to each FIS variable
    fis = AddMfToFis(fis, mf_parameters,mf_types_cell,fuzzyVarParams,fuzzyVarNames,mf_names); 
   
    % add the rule-base of the FIS
    rule_base = genRuleBase(fuzzyVarParams,rules);
    fis = addrule(fis, rule_base); 
    out = fis; 
end
