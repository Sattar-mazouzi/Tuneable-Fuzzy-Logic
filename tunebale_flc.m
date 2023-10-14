function out = tunebale_flc(input1,input2,output,mf_type ) 

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
mf_types = mf_type;
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