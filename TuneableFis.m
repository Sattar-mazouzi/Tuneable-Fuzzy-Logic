classdef  TuneableFis
    properties
     input1 ;  
     input2; 
     output;
     Rule_Base;
     Mf_Types; 
     Fis_Parameters;
     Fis;  
    end 
    methods
        function obj = TuneableFis(fis_in1, fis_in2, fis_out,fis_rules)
            if nargin > 0 
                obj.input1 = fis_in1; 
                obj.input2 = fis_in2; 
                obj.output = fis_out; 
                obj.Rule_Base = fis_rules; 
                obj.Mf_Types =  {[obj.input1.MfTypes], [obj.input2.MfTypes], [obj.output.MfTypes]};
                obj.Fis_Parameters = [obj.input1.range, obj.input1.MfNumber;   
                obj.input2.range, obj.input2.MfNumber;
                obj.output.range, obj.output.MfNumber];
            else
                disp("you must specify the Fis parameters"); 
            end
        end
        function obj = create_fis(obj)
            rules = obj.Rule_Base; 
            mftypes = obj.Mf_Types; 
            fuzzyVarParams = obj.Fis_Parameters; 
            mf_parameters   = obj.genUniformMfs(fuzzyVarParams,mftypes); 
            fis = tunebale_flc(obj.input1, obj.input2, obj.output,mftypes, rules,mf_parameters) ;
            obj.Fis = fis;
        end
        function fis = TuneMfTypes(obj)
            % Fuzzy variable paramters 
            fuzzyVarParams = obj.Fis_Parameters;  
        if nargin  >=1
            [~, mf_types_cell,~] = obj.gen_data(obj,obj.input1, obj.input2, obj.output);

            mf_params   = obj.genUniformMfs(fuzzyVarParams,mf_types_cell);
            fis = tunebale_flc(obj.input1, obj.input2, obj.output,mf_types_cell, obj.Rule_Base,mf_params) ;
    
        end 
        end 
        function output = TuneMfs(obj)  

                fuzzyVarParams = obj.Fis_Parameters;
                [mf_params ,mf_types,~] = obj.gen_data(obj,obj.input1, obj.input2, obj.output);
                mf_params  = obj.transformToMfs(fuzzyVarParams, mf_params,mf_types);
                output = tunebale_flc(obj.input1, obj.input2, obj.output,mf_types, obj.Rule_Base,mf_params) ;
        
          
        end
        function showfis(obj)
            fuzzy(obj.create_fis().Fis);
        end 
        function [mf_params, mf_types_cell, rules_data] = gen_data(obj, in1, in2, out)

            %% genrate test optimization data 
            code_bin = 3; 
            rule_row_nums = in1.MfNumber*in2.MfNumber; 
            mf_nums = in1.MfNumber+in2.MfNumber+out.MfNumber;
        
            %generate a random mf types and a random rule matrix
            mf_type_rules = randi([0,1], 1,mf_nums+(code_bin+1)*rule_row_nums);   
            mf_types = mf_type_rules(1:mf_nums);
            in1_mf_num  = in1.MfNumber; 
            in2_mf_num = in2.MfNumber;
           % out_mf_num = out.MfNumber;
            mfTypes_cell = {[mf_types(1:in1_mf_num)],[mf_types(in1_mf_num+1:in1_mf_num+in2_mf_num)], ...
            [mf_types(in1_mf_num+in2_mf_num+1:mf_nums)]} ;
            rule_data = mf_type_rules(mf_nums+1:end);
        
            % generate a random MFs params  
            trapmf_num = sum(mf_types);
            triangmf_num = sum(mf_types ==0); 
            in1Mf  = mf_types(1:in1.MfNumber); 
            in2Mf = mf_types(in1.MfNumber+1:in1.MfNumber+in2.MfNumber);
            outMf = mf_types(in1.MfNumber + in2.MfNumber+1:mf_nums); 
            in1Mf_params = 4*sum(in1Mf) + 3*sum(in1Mf==0);
            in2Mf_params = 4*sum(in2Mf) + 3*sum(in2Mf==0);
            outMf_params = 4*sum(outMf) + 3*sum(outMf==0);
        
            nVar = 3*triangmf_num + 4*trapmf_num  ; % Number Of Unkown Variable 
            VarMin = 0; % Lower Boun of the Variabls 
            VarMax =  [in1.range(2)*ones(1,in1Mf_params),...
                   in2.range(2)*ones(1,in2Mf_params),  ...
                   out.range(2)*ones(1,outMf_params)]; % Upper Bound of Decision Variable 
            VarSize = [1 nVar]; 
            rndm_mf = unifrnd(VarMin, VarMax, VarSize);
            rndm_mf_cell = {[rndm_mf(1:in1Mf_params)],...
                    [rndm_mf(in1Mf_params+1:in1Mf_params+in2Mf_params)], ... 
                    [rndm_mf(in1Mf_params+in2Mf_params+1:nVar)]};
        
            mf_params = rndm_mf_cell; 
            mf_types_cell = mfTypes_cell; 
            rules_data = rule_data; 
        end
        % this function transform a random paramaters to an accaptable Mf Parameters
        function out = transformToMfs(fisVarParameters, mf_parameters, mf_types)
            mfParams_adjusted = mf_parameters; 
            size_fisvars = size(fisVarParameters); 
            var_nums = size_fisvars(1);  % The number of Fis variables 
            
            
            for fis_var = 1:var_nums
                % Extract the parameter of MFs for the fuzzy variable
                var_range_upper = fisVarParameters(fis_var,2); 
                var_range_lower = fisVarParameters(fis_var,1); 
                mf_nums = fisVarParameters(fis_var,3);
                base =  var_range_upper/(mf_nums*4);
                mf_centers = zeros(1,mf_nums); 
        
                var_mfParams = mf_parameters{fis_var}; 
                var_mfTypes = mf_types{fis_var}; 
                idx =1; 
                for mf=1:mf_nums  % iterate through the number of mfs
                     upper_bd = var_range_upper - 2*base*(mf_nums-mf); % th upper limit of the mf 
                     lower_bd = 2*mf*base; 
                    if var_mfTypes(mf)== 0 
                        tparams = var_mfParams(idx:idx+2);
                        if mf ==1  % chechk if it is a Triangular mf
                             b = var_range_lower;   
                             a = var_range_lower; 
                             c = max(min(tparams(3),upper_bd),b+base); % set the upper bound and lowe bound of the Mf
                        elseif mf == mf_nums
                            b = var_range_upper; 
                            c = var_range_upper; 
                            a = max(min(tparams(1),var_range_upper-base),lower_bd);
                        else
                            b = max(min(tparams(2),upper_bd),mf_centers(mf-1)+base); % set the upper limit and the lower limmits;  
                            a = max(min(tparams(1),b-base),lower_bd);
                            c = max(min(tparams(3),upper_bd),b+base);
                        end
                        mf_centers(mf) = b; 
                        mfParams_adjusted{fis_var}(idx:idx+2) = [a b c]; 
                        idx = idx+3;
                    elseif var_mfTypes(mf) ==1  % if it is a trapezodial mf
                        tparams = var_mfParams(idx:idx+3);
                        if mf ==1
                            b = var_range_lower; 
                            a = var_range_lower; 
                            c = max(min(tparams(3),upper_bd),b+base); 
                            d = max(min(tparams(4),upper_bd),c+base); 
        
                        elseif mf == mf_nums
                            c = var_range_upper;  % !!
                            d = var_range_upper;  % 1!
                            b = max(min(tparams(2),d-base),mf_centers(mf-1)+base); 
                            a = max(min(tparams(1),b-base),lower_bd); %!! 
                        else
                            b = max(min(tparams(2),upper_bd),mf_centers(mf-1)+base); 
                            c = max(min(tparams(3),upper_bd),b+base); 
                            a = max(min(tparams(1),b-base),lower_bd); 
                            d = max(min(tparams(4), upper_bd),c+base);
                        end
                           mf_centers(mf) = c;
                           mfParams_adjusted{fis_var}(idx:idx+3) = [a b c d]; 
                           idx = idx+4;
                    end
                end
            end
           out = mfParams_adjusted; 
        end 
        % this function to gemerate a uniformly distributed  memebership functions 
        function out1 = genUniformMfs(obj,fuzzyVarParams, mf_types)
            mfs = {[],[],[]}; 
            mf_types_cell = mf_types; 
            varit = size(fuzzyVarParams);               % get the number of fuzzy variable 
            
            for i = 1:varit(1)                          % iterate through the number of Variables
                idx = 1;
                fuzzyVar = fuzzyVarParams(i,:); % get variable parameter
                step = fuzzyVar(2)/(fuzzyVar(3)-1);
                mf_type = mf_types_cell{i}; 
                for j = 1:fuzzyVar(3)                    % iterate through the number of MFs 
                    indic = step*j; 
                    if mf_type(j) ==0
                        mfparams = [indic-2*step, indic-step, indic];
                        mfs{i}(idx:idx+2) = mfparams ; 
                        idx = idx +3; 
                    elseif mf_type(j) == 1 
                        mfparams = [indic-2*step, indic-1.4*step,indic-0.6*step, indic];
                        mfs{i}(idx:idx+3) = mfparams ; 
                        idx = idx +4;
                    end
            
                end
            
            end
            out1 = mfs ;
            end 
        % this function generate the names of fuzzy subsets of all variables 
        function out = genMfNames(obj,mftypes) 
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
       % this function generats a a rule base for the fis 
       function out = genRuleBase(obj,fisVarParameters,rules_data)

        % get the number of memebership function for each fuccu variable
        in1_mf = fisVarParameters(1,3) ;
        in2_mf = fisVarParameters(2,3) ;
        out_mf = fisVarParameters(3,3) ;
        rules_row_nums = in1_mf*in2_mf;
        % intialize an empty rule-base matrix
        rule_base = zeros(in1_mf*in2_mf,5); 
        
        % Initialize an empty matrix to store combinations 
        rules_combinations = zeros(rules_row_nums,2);
        
        % Generate combinations using nested loops
        idx = 1; 
        for i = 1:in1_mf
            for j = 1:in2_mf
                rules_combinations(idx,:) = [i, j]; 
            idx = idx +1;
            end
        end
         
        % get the optimization data 
        output_rules = rules_data(1:3*rules_row_nums);
        %reshape the binary 1-D vector into a 2-D matrix of binary elements
        bin_rules = reshape(output_rules,rules_row_nums,3); 
        
        opperators= rules_data((3*rules_row_nums+1):end); %get the operator vector 
        %convert operator vector into ones and twos for END/OR operators
        opperators(opperators==0) = 2; 
        
        % convert the binary rules to a decimal rules 
        rules_to_dec = binaryVectorToDecimal(bin_rules); 
        % Apply rule constraints 
        %vec1_max = [ 3; 3; 3; 5; 5; 3; 3; 3; 5; 5; 3; 3; 3; 5; 5; 3; 3; 3; 5; 5; 3; 3; 3; 4; 4]; 
        %vec1_min = [ 1; 1; 1; 2; 2; 1; 1; 1; 2; 2; 1; 1; 1; 2; 2; 1; 1; 1; 2; 2; 1; 1; 1; 2; 2];
        rules_to_dec  = min(rules_to_dec,out_mf);
        rules_to_dec  = max(rules_to_dec,1); 
        
        % Construct the rule-base
        rule_base(:,1:2) = rules_combinations; 
        rule_base(:,3) = rules_to_dec;
        rule_base(:,4) = ones(1,rules_row_nums); 
        rule_base(:,5) = opperators'; 
        
        out = rule_base; 
        end
        % this functionadd the membership functions to the FIS 
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
    end 
end