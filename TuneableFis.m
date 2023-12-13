classdef  TuneableFis < handle
    properties
        input1 ;
        input2;
        output;
        Rule_Base;
        Mf_Types;
        Fis_Parameters;
        Mfs_Parameters;
        Fis;
        Problem;
        Optimization_Data;
        Tune_Flag ;
    end
    methods
        % Use the constructor to initilaze a FIS with a gevin parameters
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
                obj.create_fis();
            else
                disp("you must specify the Fis parameters");
            end
        end
        function fis = create_fis(obj)  % this function creats a fis giving the initial parameters
            fuzzyVarParams = obj.Fis_Parameters;
            mf_parameters   = obj.genUniformMfs(fuzzyVarParams,obj.Mf_Types);
            obj.Mfs_Parameters = mf_parameters;
            obj.Fis = tunebale_flc(obj.input1, obj.input2, obj.output,obj.Mf_Types, obj.Rule_Base,obj.Mfs_Parameters) ;
            fis = obj.Fis;
        end
        function output = TuneMfTypes(obj, problem) % this function tunes the types of MFs to be either a Trap mf or Trian mf
            obj.Tune_Flag ="TuneMfTypes";  %set the Tune flag to type
            disp("----> Start MF types optimization <----");
            % defging the number of decigion variable and its size and its upperand lower limits
            decvar.nvar = obj.input1.MfNumber+obj.input2.MfNumber+obj.output.MfNumber;
            decvar.varmin = 0;
            decvar.varmax = 1;
            % initialize the problem parameters
            obj.Problemdef(problem, decvar);
            tic
            obj.Optimization_Data = BGa(obj.Problem, obj.Problem);
            time = toc;
            obj.Optimization_Data.timeTaken = time/60;
            obj.runoptimization(obj.Optimization_Data.bestGenome.Position);
            output = obj.Optimization_Data;
            
        end
        function output = TuneMfs(obj,problem)  % function to tune the memebership parameters
            obj.Tune_Flag = "TuneMfs" ;
            disp("----> Start MF parameters optimization <----");
            mf_types = [obj.Mf_Types{1},obj.Mf_Types{2},obj.Mf_Types{3}];
            mf_nums = obj.input1.MfNumber+obj.input2.MfNumber+obj.output.MfNumber;
            trapmf_num = sum(mf_types);
            triangmf_num = sum(mf_types ==0);
            in1Mf  = mf_types(1:obj.input1.MfNumber);
            in2Mf = mf_types(obj.input1.MfNumber+1:obj.input1.MfNumber+obj.input2.MfNumber);
            outMf = mf_types(obj.input1.MfNumber + obj.input2.MfNumber+1:mf_nums);
            in1Mf_params = 4*sum(in1Mf) + 3*sum(in1Mf==0);
            in2Mf_params = 4*sum(in2Mf) + 3*sum(in2Mf==0);
            outMf_params = 4*sum(outMf) + 3*sum(outMf==0);
            % defining the number of decigion variable and its size and its upperand lower limits
            decvar.nvar = 3*triangmf_num + 4*trapmf_num  ; % Number Of Unkown Variable
            decvar.varmin = 0; % Lower Boun of the Variabls
            decvar.varmax =  [obj.input1.range(2)*ones(1,in1Mf_params),...
                obj.input2.range(2)*ones(1,in2Mf_params),  ...
                obj.output.range(2)*ones(1,outMf_params)];
            % initialize the problem parameters
            obj.Problemdef(problem, decvar)
            tic;
            obj.Optimization_Data = RGa(obj.Problem, obj.Problem);
            time = toc;
            obj.Optimization_Data.timeTaken = time/60;
            obj.runoptimization(obj.Optimization_Data.bestGenome.Position);
            output = obj.Optimization_Data;
            
            
        end
        function output = TuneRules(obj,problem)
            disp("----> Start Rule-base optimization <----");
            obj.Tune_Flag = "TuneRules" ;
            [decvar.nvar, decvar.varmin, decvar.varmax  ] = obj.rule_data();
            obj.Problemdef(problem,decvar);
            tic;
            obj.Optimization_Data = BGa(obj.Problem, obj.Problem);
            time = toc;
            obj.Optimization_Data.timeTaken = time/60;
            obj.runoptimization(obj.Optimization_Data.bestGenome.Position);
            output = obj.Optimization_Data;
            
        end
        % function to tune the fuzzy varieble ranges
        function output = TuneFuzzyVarsRange(obj, problem)
            disp("----> Start Fuzzy Variables Ranges optimization <----");
            [decvar.nvar, decvar.varmin, decvar.varmax  ]  = obj.FVR_data();
            obj.Tune_Flag = "TuneFVR";
            obj.Problemdef(problem, decvar)
            tic;
            obj.Optimization_Data = RGa(obj.Problem, obj.Problem);
            time = toc;
            obj.Optimization_Data.timeTaken = time/60;
            obj.runoptimization(obj.Optimization_Data.bestGenome.Position);
            output = obj.Optimization_Data;
            
        end
        
        function output = TuneMfsNumbers(obj, problem)
            disp("----> Start number of MFs optimization <----");
            obj.Tune_Flag = "TuneMfNumber" ;
            [decvar.nvar, decvar.varmin, decvar.varmax  ] = obj.Mfs_nums_data();
            obj.Problemdef(problem,decvar);
            tic;
            obj.Optimization_Data = BGa(obj.Problem, obj.Problem);
            time = toc;
            obj.Optimization_Data.timeTaken = time/60;
            obj.runoptimization(obj.Optimization_Data.bestGenome.Position);
            output = obj.Optimization_Data;
          
        end
        function output = TuneAll(obj,problem)
            disp("----> Start full optimization")
            obj.TuneMfTypes(problem);
            obj.TuneRules(problem);
            output = obj.TuneMfs(problem)  ;
        end
        function runoptimization(obj,X)
            if obj.Tune_Flag == "TuneMfTypes"
                
                obj.Mf_Types = obj.mftype_data(X);
                obj.Mfs_Parameters   = obj.genUniformMfs(obj.Fis_Parameters, obj.Mf_Types );
                obj.Fis = tunebale_flc(obj.input1, obj.input2, obj.output, obj.Mf_Types, obj.Rule_Base, obj.Mfs_Parameters );
                
            elseif obj.Tune_Flag == "TuneMfs"
                
                mf_params_data = obj.mf_parameter_data(X);
                obj.Mfs_Parameters  = obj.transformToMfs(obj.Fis_Parameters, mf_params_data, obj.Mf_Types);
                obj.Fis = tunebale_flc(obj.input1, obj.input2, obj.output, obj.Mf_Types, obj.Rule_Base, obj.Mfs_Parameters) ;
                
            elseif obj.Tune_Flag == "TuneRules"
                
                obj.Rule_Base = genRuleBase(obj,X);
                obj.Fis =  tunebale_flc(obj.input1, obj.input2, obj.output,obj.Mf_Types, obj.Rule_Base, obj.Mfs_Parameters);
                
            elseif obj.Tune_Flag == "TuneFVR"
                obj.input1.range = [0 X(1)];
                obj.input2.range = [0 X(2)];
                obj.output.range = [0 X(3)];
                obj.Fis_Parameters = [obj.input1.range, obj.input1.MfNumber;
                    obj.input2.range, obj.input2.MfNumber;
                    obj.output.range, obj.output.MfNumber];
                fuzzyVarParams = obj.Fis_Parameters;
                obj.Mfs_Parameters   = obj.genUniformMfs(fuzzyVarParams,obj.Mf_Types);
                obj.Fis =   tunebale_flc(obj.input1, obj.input2, obj.output,obj.Mf_Types, obj.Rule_Base, obj.Mfs_Parameters);
            elseif obj.Tune_Flag == "TuneMfNumber"
                update_mf_numbers(obj,X);
                code_bin = 3;
                rule_row_nums = obj.input1.MfNumber*obj.input2.MfNumber;
                nvar_rule = (code_bin+1)*rule_row_nums ;
                random_rules  = randi([0 1],[1 nvar_rule]);
                obj.Rule_Base = genRuleBase(obj,random_rules);
                mf_parameters   = obj.genUniformMfs(obj.Fis_Parameters,obj.Mf_Types);
                obj.Mfs_Parameters = mf_parameters;
                obj.Fis =  tunebale_flc(obj.input1, obj.input2, obj.output,obj.Mf_Types, obj.Rule_Base, obj.Mfs_Parameters);
                
                
                
                
            end
        end
        function showfis(obj)
            fuzzy(obj.Fis);
        end
        % this function generate a random data
        function [mf_params, mf_types_cell, rules_data] = gen_data(obj,in1, in2, out)
            
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
            
            % generate a random MFs obj.Problem
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
        % This function to calclate the length of 1-d array of rules
        function [nvar, varmin, varmax]= rule_data(obj)
            %% genrate test optimization data
            code_bin = 3;
            rule_row_nums = obj.input1.MfNumber*obj.input2.MfNumber;
            nvar = (code_bin+1)*rule_row_nums ;
            varmin = 0;
            varmax = 1;
        end
        % this function generats the data for optimizaing the fuzzy variable ranages
        function [nvar, varmin, varmax] = FVR_data(obj)
            nvar = 3;
            varmin = [obj.input1.range(2)/3, obj.input2.range(2)/3, obj.output.range(2)/3];
            varmax = [obj.input1.range(2), 2*obj.input2.range(2), 2*obj.output.range(2)] ;
            
        end
        % this function generates the data for tuning the Mfs Number tuning
        function [nvar, varmin, varmax] = Mfs_nums_data(obj)
            %MaxMfsNumber = 12;
            bit_number  = 4;  % number of bits to code the decision variabel range
            %fuzzy_var_nums = 3; % number of fuzzy variables
            %dec_minvar= 3;      % minimum number of memebership function for all variables
            %binary_minvar = decimalToBinaryVector(dec_minvar,bit_number);
            %binary_maxvar = decimalToBinaryVector(MaxMfsNumber,bit_number);
            
            nvar = 3*bit_number ;
            %varmin = repmat(binary_minvar,1, fuzzy_var_nums) ;
            %varmax =  repmat(binary_maxvar, 1, fuzzy_var_nums);
            varmin = 0;
            varmax = 1;
        end
        function  update_mf_numbers(obj, X)
            reshaped_X = reshape(X,4,3)';
            mf_numbers = binaryVectorToDecimal(reshaped_X)';
            % update the mf numbers of each  fuzzy variable
            mf_numbers  = min(mf_numbers,12);
            mf_numbers  = max(mf_numbers,3);
            obj.input1.MfNumber = mf_numbers(1);
            obj.input2.MfNumber = mf_numbers(2);
            obj.output.MfNumber = mf_numbers(3);
            
            % generate a random mf types corresponding to the new mf numbers
            obj.input1.MfTypes = randi([0 1],[1 obj.input1.MfNumber]);
            obj.input2.MfTypes = randi([0 1],[1 obj.input2.MfNumber]);
            obj.output.MfTypes = randi([0 1],[1 obj.output.MfNumber]);
            
            % update class params
            obj.Mf_Types =  {[obj.input1.MfTypes], [obj.input2.MfTypes], [obj.output.MfTypes]};
            obj.Fis_Parameters = [obj.input1.range, obj.input1.MfNumber;
                obj.input2.range, obj.input2.MfNumber;
                obj.output.range, obj.output.MfNumber];
            
        end
        
        function output = mftype_data(obj, mf_types)
            
            mf_nums = obj.input1.MfNumber+obj.input2.MfNumber+obj.output.MfNumber;
            in1_mf_num  = obj.input1.MfNumber;
            in2_mf_num = obj.input2.MfNumber;
            % out_mf_num = out.MfNumber;
            mfTypes_cell = {[mf_types(1:in1_mf_num)],[mf_types(in1_mf_num+1:in1_mf_num+in2_mf_num)], ...
                [mf_types(in1_mf_num+in2_mf_num+1:mf_nums)]} ;
            output = mfTypes_cell;
        end
        % this function create a cell array of mf parameter
        function output = mf_parameter_data(obj, X)
            mf_types = [obj.Mf_Types{1},obj.Mf_Types{2},obj.Mf_Types{3}];
            % generate a random MFs obj.Problem
            mf_nums = obj.input1.MfNumber+obj.input2.MfNumber+obj.output.MfNumber;
            trapmf_num = sum(mf_types);
            triangmf_num = sum(mf_types ==0);
            in1Mf  = mf_types(1:obj.input1.MfNumber);
            in2Mf = mf_types(obj.input1.MfNumber+1:obj.input1.MfNumber+obj.input2.MfNumber);
            outMf = mf_types(obj.input1.MfNumber + obj.input2.MfNumber+1:mf_nums);
            in1Mf_params = 4*sum(in1Mf) + 3*sum(in1Mf==0);
            in2Mf_params = 4*sum(in2Mf) + 3*sum(in2Mf==0);
            % outMf_params = 4*sum(outMf) + 3*sum(outMf==0);
            
            nVar = 3*triangmf_num + 4*trapmf_num  ; % Number Of Unkown Variable
            
            rndm_mf = X;
            rndm_mf_cell = {[rndm_mf(1:in1Mf_params)],...
                [rndm_mf(in1Mf_params+1:in1Mf_params+in2Mf_params)], ...
                [rndm_mf(in1Mf_params+in2Mf_params+1:nVar)]};
            output = rndm_mf_cell;
            
        end
        
        % this function transform a random paramaters to an accaptable Mf Parameters
        function out = transformToMfs(obj,fisVarParameters, mf_parameters, mf_types)
            mfParams_adjusted = mf_parameters;
            size_fisvars = size(fisVarParameters);
            var_nums = size_fisvars(1);  % The number of Fis variables
            
            for fis_var = 1:var_nums
                % Extract the parameter of MFs for the fuzzy variable
                var_range_upper = fisVarParameters(fis_var,2);
                var_range_lower = fisVarParameters(fis_var,1);
                mf_nums = fisVarParameters(fis_var,3);
                base =  var_range_upper/(mf_nums*3);
                mf_centers = zeros(1,mf_nums);
                
                var_mfParams = mf_parameters{fis_var};
                var_mfTypes = mf_types{fis_var};
                idx =1;
                for mf=1:mf_nums  % iterate through the number of mfs
                    upper_bd = var_range_upper - 2*base*(mf_nums-mf); % th upper limit of the mf
                    lower_bd = base;
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
        % this function to generats a uniformly distributed  memebership functions
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
        % this function generats a a rule base for the fis based on the given rule date
        function out = genRuleBase(obj,rules_data)
            fisVarParameters = obj.Fis_Parameters;
            % get the number of memebership function for each fuzzy variable
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
        % this function set the parameters of the optimization function
        function  Problemdef(obj,problem,decvar)
            %% GA Parameters
            obj.Problem.CostFunction =  problem.CostFunction;
            
            select_func = ["Random", "RuletteWheel", "Hybrid"];
            crossover_func = ["SinglePoint", "DoublePoint", "Uniform", "Hybrid"];
            obj.Problem.MaxIt= problem.MaxIt;
            obj.Problem.nPop = problem.nPop;   % Population  Genome
            obj.Problem.pc = 1;      % Number of Offsprings
            obj.Problem.mu = 0.3;   % Mutation Rate
            obj.Problem.betea = 1;   % Beta value for probablilty distributions
            obj.Problem.sigma = 0.25;  % step size
            obj.Problem.gamma = 0.12;  %  cross over
            obj.Problem.SelectionFunc = select_func(2);
            obj.Problem.CrossoverFunc = crossover_func(3);
            % Decision Variable parameters
            obj.Problem.nVar = decvar.nvar;
            obj.Problem.VarMin = decvar.varmin;
            obj.Problem.VarMax = decvar.varmax;
            
        end
        
    end
end