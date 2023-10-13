function out = tunebale_flc(input1,input2,output ) 

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

fuzzyVarParams = [in1_range(1), in1_mfNumber;   % Fuzzy variable paramters 
                  in2_range(1), in2_mfNumber;
                  out_range(1), out_mfNumber;] ; 
fuzzyVarNames = [in1_name, in2_name, out_name];  % FIS variable names  

%% generate a uniformly distributed memebership function 
 
mf_parameters, mf_types_cell = genUniformMfs(mf_types, fuzzyVarParams); 






end