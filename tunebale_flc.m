function out = tunebale_flc(input1,input2,output,optimization_data ) 
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
    
    fuzzyVarParams = [in1_range(2), in1_mfNumber;   % Fuzzy variable paramters 
                      in2_range(2), in2_mfNumber;
                      out_range(2), out_mfNumber;] ; 
    fuzzyVarNames = [in1_name, in2_name, out_name];  % FIS variable names 
    mf_nums = in1_mfNumber + in2_mfNumber +  out_mfNumber; % number of all MFs
    % optimization information;
  
    mf_types = optimization_data(1:mf_nums);
    rule_data = optimization_data(mf_nums+1:end);
    %% generate a uniformly distributed memebership function
    % generate MFs parameters and MF types
    [mf_parameters, mf_types_cell] = genUniformMfs(fuzzyVarParams,mf_types); 
    
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
    rule_base = genRuleBase(fuzzyVarParams,rule_data);
    fis = addrule(fis, rule_base); 
    out = fis; 
end
