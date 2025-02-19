classdef  TuneableFis < handle
    % this Class used to tune a mamdani fuzzy logic inference system (FIS)
    % the class create a tunebale FIS, to tune a set of its parameters;
    % ARGUMENTS:
    %  FIS_IN1, FIS_IN2, FIS_IN3 :  structure-type to  discribe the fuzzy variable parametes,
    %    which contain three fields:
    %       name    : fuzzy variabel name (string)
    %       range   : a tow element array contains the upper and lowe limits of fuzzy variable
    %       MfNumber: the number of memebership functions
    %       Mftypes : 1-D binary array contains the type of each MF,0 for Triangular,1 for Trapezodial
    % the followeing parameter can be tuned:
    % 1- Fuzzy memebership function types are tuned by altering between Triangular\Trapezodial mf types
    %   using TuneMfTypes(problem) class function
    % 2- Fuzzy variable ranges are tuned using TuneFuzzyVarsRange(problem) class method
    % 3- Membership function parameters are tune using the TuneMfParameters(problem)
    % 4- The Rule base are tuned using the TuneRules(problem, arg1,arg2) method where:
    %       arg2: is a boolean variable to decide weather tuning the AND\OR operator True for tuning,False otherwise.
    %       arg3: is a cell array {[upper bounds], [Lower bounds]}
    % 5- The Rule Wights are optimized using TuneRuleWeigths() method
    % 6- The number of membership function tuned using the TuneMfsNumbers()
    %       !!!!! stil working on this, when changing the nuber of MfS the rules has to be updated, and the rules usualy are designe
    %       with human knowldge of and informations about the target system,
    %       the function tunes the memebership function number, and every time generates a random set of rules;
    %
    properties
        input1 struct;
        input2 struct;
        output struct;
        Rule_Base (:,5) double;
        Mf_Types cell;
        Fis_Parameters double;
        Mfs_Parameters cell;
        Fis;
        Problem =  struct('CostFunction',[],'MaxIt',[],'nPop',[]);
        Optimization_Data  struct;
        Tune_Flag string = "None";
        Tune_options = struct('vartype',[],'tune_operators',[],'rule_limits',[], 'tune_range',[], 'range_offset',[]) ;
        
    end
    methods
        % Use the constructor to initilaze a FIS with a gevin parameters
        function obj = TuneableFis(fis_in1, fis_in2, fis_out,fis_rules)
            arguments
                fis_in1 struct;
                fis_in2 struct;
                fis_out struct;
                fis_rules  (:,5) double ;
            end
            if nargin > 0
                obj.input1 = fis_in1;
                obj.input2 = fis_in2;
                obj.output = fis_out;
                obj.Rule_Base = fis_rules;
                obj.Mf_Types =  {[obj.input1.MfTypes], [obj.input2.MfTypes], [obj.output.MfTypes]};
                obj.Fis_Parameters = [obj.input1.range, obj.input1.MfNumber;
                    obj.input2.range, obj.input2.MfNumber;
                    obj.output.range, obj.output.MfNumber];
                obj.Tune_options.tune_operators  = false;
                %obj.Tune_Flag = "None" ;
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
            if not(obj.Tune_Flag == "TuneMfs" )
                obj.Tune_Flag ="TuneMfTypes";  %set the Tune flag
                obj.Tune_options.vartype = 'bitstring';
                disp("----> Start MF types optimization ");
                % define the number of decision variable and its size and its upperlimit lowerlimits
                decvar.nvar = obj.input1.MfNumber+obj.input2.MfNumber+obj.output.MfNumber;
                decvar.varmin = [];
                decvar.varmax = [];
                % initialize the problem parameters
                obj.Optimization_Data  = start_optimize(obj, problem, decvar );
                obj.evaluate(obj.Optimization_Data.sol);
                output = obj.Optimization_Data;
            else
                disp("Enable to tune the membership function types after tuning membership function parameters.")
            end
            
        end
        function output = TuneMfParameters(obj,problem)  % function to tune the memebership parameters
            obj.Tune_Flag = "TuneMfs" ;
            disp("----> Start MF parameters optimization ");
            obj.Tune_options.vartype = 'doubleVector';
            
            [decvar.nvar, decvar.varmin, decvar.varmax]=  obj.mfParamsData();
            % initialize the problem parameters
            obj.Optimization_Data = start_optimize(obj, problem, decvar );
            obj.evaluate(obj.Optimization_Data.sol);
            
            output = obj.Optimization_Data;
            
            
        end
        function output = TuneRules(obj,problem, tune_operator, rule_limits)
            arguments
                obj TuneableFis;
                problem struct;
                tune_operator logical = false;
                rule_limits  double = [obj.output.MfNumber; 1];
            end
            obj.Tune_options.tune_operators = tune_operator;
            obj.Tune_options.vartype = 'bitstring';
            
            obj.Tune_options.rule_limits = rule_limits;
            if not (obj.Tune_Flag == "TuneRW")
                disp("----> Start Rule-base optimization <----");
                obj.Tune_Flag = "TuneRules" ;
                [decvar.nvar, decvar.varmin, decvar.varmax  ] = obj.rule_data();
                decvar.varmin = [] ;
                decvar.varmax = [];
                obj.Problemdef(problem,decvar);
                obj.Optimization_Data  = start_optimize(obj, problem, decvar );
                obj.evaluate(obj.Optimization_Data.sol);
                output = obj.Optimization_Data;
            else
                disp("Rule must be tuned before tuning Rule weights")
            end
            
        end
        %this fucntion tunes the rules wieghts
        function output = TuneRuleWeigths(obj, problem)
            disp("----> Start ruel weights optimization ");
            obj.Tune_Flag = "TuneRW";
            obj.Tune_options.vartype = 'doubleVector';
            [nvar, varmin, varmax] = obj.ruleWeightsdata();
            decvar.nvar = nvar;
            decvar.varmin = varmin;
            decvar.varmax =  varmax;
            % obj.Problemdef(problem, decvar)
            %    tic;
            %    obj.Optimization_Data = RGa(obj.Problem, obj.Problem);
            %    time = toc;
            %    obj.Optimization_Data.timeTaken = time/60;
            %    obj.evaluate(obj.Optimization_Data.bestGenome.Position);
            obj.Optimization_Data  = start_optimize(obj, problem, decvar );
            obj.evaluate(obj.Optimization_Data.sol);
            
            output = obj.Optimization_Data;
            
        end
        % function to tune the fuzzy varieble ranges
        function output = TuneFuzzyVarsRange(obj, problem, tune_range,range_offset)
            arguments
                obj TuneableFis;
                problem struct;
                tune_range string = "upper";
                range_offset double = 0.2;
            end
            obj.Tune_options.tune_range = tune_range;
            obj.Tune_options.range_offset =range_offset;
            obj.Tune_options.vartype = 'doubleVector';
            
            if not(obj.Tune_Flag == "TuneMfs")
                disp("----> Start Fuzzy Variables Ranges optimization ");
                [decvar.nvar, decvar.varmin, decvar.varmax  ]  = obj.FVR_data();
                obj.Tune_Flag = "TuneFVR";
                obj.Problemdef(problem, decvar)
                %tic;
                %obj.Optimization_Data = RGa(obj.Problem, obj.Problem);
                %time = toc;
                %obj.Optimization_Data.timeTaken = time/60;
                %
                % obj.evaluate(obj.Optimization_Data.bestGenome.Position);
                obj.Optimization_Data  = start_optimize(obj, problem, decvar );
                obj.evaluate(obj.Optimization_Data.sol);
                output = obj.Optimization_Data;
            else
                disp("Enable to tune the fuzzy variablee range ofter tuning the MFs parameters")
            end
            
        end
        
        function output = TuneMfsNumbers(obj, problem)
            if (obj.Tune_Flag == "None")
                disp("----> Start number of MFs optimization ");
                obj.Tune_Flag = "TuneMfNumber" ;
                [decvar.nvar, decvar.varmin, decvar.varmax  ] = obj.Mfs_nums_data();
                obj.Problemdef(problem,decvar);
                tic;
                obj.Optimization_Data = BGa(obj.Problem, obj.Problem);
                time = toc;
                obj.Optimization_Data.timeTaken = time/60;
                obj.evaluate(obj.Optimization_Data.bestGenome.Position);
                output = obj.Optimization_Data;
            else
                disp("Enable to tune the number of MF")
            end
        end
        function output = TuneAll(obj,problem)
            if not(obj.Tune_Flag == "TuneMfs" )
                disp("----> Start full optimization")
                it = 0;
                while it < 5
                    if it == 0
                        mkdir 1.Tune_Mf_types ;
                        cd 1.Tune_Mf_types ;
                        out.mftypes = obj.TuneMfTypes(problem);
                    elseif it ==1
                        mkdir 2.Tune_Fuzzy_Ranges
                        cd 2.Tune_Fuzzy_Ranges
                        out.fvrange = obj.TuneFuzzyVarsRange(problem);
                    elseif it == 2
                        mkdir 3.Tune_Mfs
                        cd 3.Tune_Mfs
                        out.mfparams = obj.TuneMfParameters(problem);
                    elseif it == 3
                        mkdir 4.Tune_Rule_Weights
                        cd 4.Tune_Rule_Weights
                        out.rulebase = obj.TuneRules(problem);
                        
                    elseif it == 4
                        mkdir 5.Tune_Rulebase
                        cd 5.Tune_Rulebase
                        out.ruleweights = obj.TuneRuleWeigths(problem);
                    end
                    cd ..;
                    it = it +1;
                end
                output = out;
                
            else
                disp("Enable to prform  the tune all after tuning membership function parameters.")
            end
        end
        function output = Optimizefis(obj, problem)
            % This function optimizes the FIS rule base, memebership functions, and rule weights
            arguments
                obj TuneableFis;
                problem struct;
            end
            
            disp("----> Start comperhensive optimization ");
            obj.Tune_Flag = "CompOpt" ;
            obj.Tune_options.vartype = 'doubleVector';
            % define the number of decision variable and its size and its upperand lower limits
            % mf parameters, rules, rule weights
            [nvar, varmin, varmax] = obj.optimizeFisData();
            decvar.nvar = nvar;
            decvar.varmin = varmin;
            decvar.varmax = varmax;
            obj.Problemdef(problem,decvar);
            
            obj.Optimization_Data  = start_optimize(obj, problem, decvar );
            obj.evaluate(obj.Optimization_Data.sol);
            disp("----> Optmization completed ");
            output = obj.Optimization_Data;
            
            
        end
        function evaluate(obj,X)
            arguments
                obj TuneableFis;
                X double;
            end
            %disp("Solution")
            %disp(X);
            [var_ranges, mf_types, mf_params, rule_base, rule_connections, rule_weights ] = get_fis_parameters(obj);
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
                obj.input1.range = [X(1) X(4)];
                obj.input2.range = [X(2) X(5)];
                obj.output.range = [X(3) X(6)];
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
            elseif obj.Tune_Flag == "TuneRW"
                rule_weights = X' ;
                obj.Rule_Base(:,4) = rule_weights;
                obj.Fis =  tunebale_flc(obj.input1, obj.input2, obj.output,obj.Mf_Types, obj.Rule_Base, obj.Mfs_Parameters);
            elseif obj.Tune_Flag == "CompOpt"
                % update the FIS parameters
                size_rule = size(rule_base(:,3)');
                size_mfparams = size(mf_params);
                size_rw = size(rule_weights);
                x_mf = X(1:size_mfparams(2));
                x_rule = X(size_mfparams(2)+1:size_mfparams(2)+size_rule(2));
                x_rule_weights  = X(size_mfparams(2)+size_rule(2)+1:end);
                
                mf_params_data = obj.mf_parameter_data(x_mf);
                obj.Mfs_Parameters  = obj.transformToMfs(obj.Fis_Parameters, mf_params_data, obj.Mf_Types);
                obj.Rule_Base(:,3)= x_rule';
                obj.Rule_Base(:,4) = x_rule_weights';
                obj.Fis =  tunebale_flc(obj.input1, obj.input2, obj.output,obj.Mf_Types, obj.Rule_Base, obj.Mfs_Parameters);
                
            end
        end
        function showfis(obj)
            fuzzy(obj.Fis);
        end
        
        function [mf_params, mf_types_cell, rules_data] = gen_data(obj,in1, in2, out)
            % This function generate a random data
            % genrate test optimization data
            % the function used for testing, when there is no real optimization
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
        function [nvar, varmin, varmax] = optimizeFisData(obj)
            %This function calculate the data for optimizefis(), it calculatet he number of decision variables, their upper and lower limits
            [nvar1, varmin1, varmax1]=  obj.mfParamsData();
            
            nvar2 = obj.input1.MfNumber*obj.input2.MfNumber;
            varmin2 = ones(1,nvar2);
            varmax2 =  obj.output.MfNumber * ones(1,nvar2);
            
            [nvar3, varmin3, varmax3] = obj.ruleWeightsdata();
            
            
            nvar = nvar1 + nvar2 + nvar3;
            varmin = [varmin1 varmin2 varmin3];
            varmax = [varmax1 varmax2 varmax3];
            
        end
        
        function [nvar, varmin, varmax]=  mfParamsData(obj)
            %This function calculate the data for TuneMfParameters(), it calculatet he number of decision variables, their upper and lower limits
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
            nvar = 3*triangmf_num + 4*trapmf_num  ; % Number Of Unkown Variable
            varmin = zeros(1,nvar); % Lower Boun of the Variabls
            varmax =  [obj.input1.range(2)*ones(1,in1Mf_params),...
                obj.input2.range(2)*ones(1,in2Mf_params),  ...
                obj.output.range(2)*ones(1,outMf_params)];
        end
        function [nvar, varmin, varmax] = ruleWeightsdata(obj)
            %This function calculate the data for TuneRuleWeigths(), it calculatet he number of decision variables, their upper and lower limits
            nvar = obj.input1.MfNumber*obj.input2.MfNumber;
            varmin = zeros(1,nvar);
            varmax =  ones(1,nvar);
            
        end
        
        
        function [nvar, varmin, varmax]= rule_data(obj)
            % This function to calclate the length of 1-d array of rules
            % for tuning the the rule-base
            %% genrate test optimization data
            code_bin = 3;
            rule_row_nums = obj.input1.MfNumber*obj.input2.MfNumber;
            nvar = (code_bin+1)*rule_row_nums ;
            varmin = 0;
            varmax = 1;
        end
        % this function generats the data for optimizaing the fuzzy variable ranages
        function [nvar, varmin, varmax] = FVR_data(obj)
            nvar = 6;
            offset_rate = obj.Tune_options.range_offset;
            %size_fisvars = size(fisVarParameters);
            %var_nums = size_fisvars(1);
            
            % get the range length of each variable
            norm_range1 = obj.normrange(obj.Fis_Parameters(1,1:2));
            var_len1 = norm_range1(2);
            norm_range2 = obj.normrange(obj.Fis_Parameters(2,1:2));
            var_len2 = norm_range2(2);
            norm_range3 = obj.normrange(obj.Fis_Parameters(3,1:2));
            var_len3 = norm_range3(2);
            %calculate the offset
            offset1 = offset_rate*var_len1;
            offset2 = offset_rate*var_len2;
            offset3 = offset_rate*var_len3;
            
            if obj.Tune_options.tune_range == "both"
                
                varmin = [(obj.input1.range(1) - offset1), (obj.input2.range(1) -offset2 ), (obj.output.range(1) - offset3), ...
                    (obj.input1.range(2) -offset1), (obj.input2.range(2) - offset2), (obj.output.range(2) -offset3)];
                
                varmax = [(obj.input1.range(1) + offset1), (obj.input2.range(1) + offset2), (obj.output.range(1) +offset3), ...
                    (obj.input1.range(2) + offset1), (obj.input2.range(2) + offset2), (obj.output.range(2) +offset3)] ;
            elseif obj.Tune_options.tune_range  == "lower"
                
                varmin = [(obj.input1.range(1) - offset1), (obj.input2.range(1) -offset2 ), (obj.output.range(1) - offset3), ...
                    (obj.input1.range(2)), (obj.input2.range(2)), (obj.output.range(2))];
                
                varmax = [(obj.input1.range(1) + offset1), (obj.input2.range(1) + offset2), (obj.output.range(1) +offset3), ...
                    (obj.input1.range(2)), (obj.input2.range(2)), (obj.output.range(2))] ;
                
            else
                varmin = [(obj.input1.range(1)), (obj.input2.range(1)), (obj.output.range(1) ), ...
                    (obj.input1.range(2) -offset1), (obj.input2.range(2) - offset2), (obj.output.range(2) -offset3)];
                
                varmax = [(obj.input1.range(1)), (obj.input2.range(1)), (obj.output.range(1)), ...
                    (obj.input1.range(2) + offset1), (obj.input2.range(2) + offset2), (obj.output.range(2) +offset3)] ;
                
            end
            
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
        % this function updates the new number of memebership functions
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
            offsets = zeros(1,3);
            for fis_var = 1:var_nums
                % Extract the parameter of MFs for the fuzzy variable
                range = fisVarParameters(fis_var,1:2);
                norm_range = obj.normrange(range);
                %disp(norm_range(1:2));
                offsets(fis_var) = norm_range(3);
                %disp(offset);
                %var_range_upper = fisVarParameters(fis_var,2);
                %var_range_lower = fisVarParameters(fis_var,1);
                var_range_upper = norm_range(2);
                var_range_lower = norm_range(1);
                
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
            
            out = obj.shiftMfparameters(mfParams_adjusted,offsets);
        end
        % this function to generats a uniformly distributed  memebership functions
        function out1 = genUniformMfs(obj,fuzzyVarParams, mf_types)
            mfs = {[],[],[]};                           % Initialize an empty cell array to store membership function parameters
            mf_types_cell = mf_types;
            varit = size(fuzzyVarParams);               % get the number of fuzzy variable
            
            for i = 1:varit(1)                          % iterate through the number of FIS Variables
                idx = 1;                                % set index to keep track for membership functions parameters
                fuzzyVar = fuzzyVarParams(i,:);         % get ech variable parameter [range, mf number]
                range = [fuzzyVar(1), fuzzyVar(2)];
                %disp(range)
                paramss = obj.normrange(range) ;
                %disp(paramss);
                norm_range  = paramss(1:2); %normalized range
                offset = paramss(3);
                % offset = params(3);
                step = norm_range(2)/(fuzzyVar(3)-1);     % range/(mfnumber -1)
                mf_type = mf_types_cell{i};
                for j = 1:fuzzyVar(3)                    % iterate through the number of MFs
                    indic = step*j;
                    if mf_type(j) ==0
                        mfparams = [indic-2*step, indic-step, indic] - offset;
                        mfs{i}(idx:idx+2) = mfparams ;
                        idx = idx +3;
                    elseif mf_type(j) == 1
                        mfparams = [indic-2*step, indic-1.4*step,indic-0.6*step, indic] - offset;
                        mfs{i}(idx:idx+3) = mfparams ;
                        idx = idx +4;
                    end
                    
                end
                
            end
            out1 = mfs ;
        end
        function out =  normrange(obj,range)
            % this function normalize range to 0 to 1
            fst = range(1) ;
            scnd = range(2);
            lw = 0;
            up = 0;
            offest = 0;
            
            if  fst <=0
                offest  = abs(fst);
                %lw = fst +abs(fst);
                %up = scnd + abs(fst);
                
            elseif fst > 0
                offest = -fst;
                %lw = fst - fst;
                %up = scnd - fst;
            end
            lw = fst + offest ;
            up = scnd + offest;
            out = [lw, up, offest] ;
            
        end
        function output = shiftMfparameters(obj, mf_parameters, offset)
            %% this function shifts the generated Mf parameter by the given offset
            size_cellarray = size(mf_parameters);
            for i =1:size_cellarray(2)              % itereate throgh the number of fuzzy variable
                mf_parameters{i} = mf_parameters{i} - offset(i);
            end
            output = mf_parameters;
        end
        % this function generats a a rule base for the fis based on the given rule date
        function out = genRuleBase(obj,rules_data)
            arguments
                obj TuneableFis;
                rules_data double
            end
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
            %get the operator vector
            opperators= rules_data((3*rules_row_nums+1):end);
            if obj.Tune_options.tune_operators  == true
                %convert operator vector into ones and twos for END/OR operators
                opperators(opperators==0) = 2;
            else
                % make the operator vector all ones
                opperators(opperators==0) = 1;
            end
            
            % convert the binary rules to a decimal rules
            rules_to_dec = binaryVectorToDecimal(bin_rules);
            
            % Apply rule constraints
            if isequal(obj.Tune_options.rule_limits, [obj.output.MfNumber; 1])
                rules_to_dec  = min(rules_to_dec,out_mf);
                rules_to_dec  = max(rules_to_dec,1);
            else
                %disp(rules_to_dec)
                rules_to_dec  = min(rules_to_dec,obj.Tune_options.rule_limits(1,:)');
                rules_to_dec  = max(rules_to_dec,obj.Tune_options.rule_limits(2,:)');
            end
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
        %%%%
        function output = start_optimize(obj, problem, decvar )
            
            obj.Problemdef(problem, decvar);
            MaxGens = obj.Problem.MaxIt  ;
            PopSize = obj.Problem.nPop;
            %disp("Varmin");
            %disp(obj.Problem.VarMin ) ;
            %disp("VarMax")
            %disp(obj.Problem.VarMax) ;
            [var_ranges, mf_types, mf_params, rule_base, rule_connections, rule_weights ] = get_fis_parameters(obj);
            initial_var = zeros(1,decvar.nvar);
            
            
            if obj.Tune_Flag == "TuneMfTypes"
                initial_var = mf_types;
                intcon = [] ;
            elseif obj.Tune_Flag == "TuneMfs"
                initial_var = mf_params;
                intcon = [] ;
            elseif obj.Tune_Flag == "TuneRules"
                rule_consequent = rule_base(:,3);
                connections  = rule_connections;
                a_bin = []; % Initialize empty array for binary digits
                for i = 1:length(rule_consequent)
                    % Convert each number to a 3-bit binary string
                    binary_str = dec2bin(rule_consequent(i), 3);
                    % Convert the binary string to numeric digits and append to the array
                    a_bin = [a_bin, str2num(binary_str(:))']; % Convert string to double and append
                end
                initial_var = [a_bin connections];
                intcon = [] ;
                
            elseif obj.Tune_Flag == "TuneFVR"
                initial_var = var_ranges;
                intcon = [] ;
            elseif obj.Tune_Flag == "TuneMfNumber"
                initial_var = zeros(1,decvar.nvar);
                intcon = [] ;
            elseif obj.Tune_Flag == "TuneRW"
                initial_var = rule_weights;
                intcon = [] ;
            elseif obj.Tune_Flag == "CompOpt"
                initial_var = [mf_params rule_base(:,3)' rule_weights];
                disp("Initial Var" + num2str(size(initial_var)));
                disp("Nvar= " + num2str(decvar.nvar));
                
                size_rule = size(rule_base(:,3)');
                size_mfparams = size(mf_params);
                startt = size_mfparams(2) + 1;
                endd = size_mfparams(2) + size_rule(2);
                intcon = [startt:1:endd];
            end
            
            % Initialize a matrix to store the initial population
            initialPop = [];
            initialPop(1,:) = initial_var;
            options = optimoptions('ga', PopulationType = obj.Tune_options.vartype ,PlotFcn = @gaplotbestf, OutputFcn = @captureBestCosts, MaxGenerations = MaxGens,PopulationSize = PopSize, InitialPopulationMatrix = initialPop);
            [x,fval,exitflag,out, population,scores] = ga(obj.Problem.CostFunction, decvar.nvar,[],[], [],[],...
                obj.Problem.VarMin,obj.Problem.VarMax, [], intcon, options);
            
            output.sol = x ;
            output.best = fval;
            output.exitflag = exitflag;
            output.out = out;
            output.population = population;
            output.scores = scores;
            
        end
        function [var_ranges, mf_types, mf_params, rule_base,rule_connections, rule_weights ] = get_fis_parameters(obj);
            %% Get Fis paramater for optimization
            fis = obj.Fis;
            inputs = fis.Inputs;
            ouputs = fis.Outputs ;
            
            %% Extract mf types and parameters;
            input_mfs = {inputs.MembershipFunctions};
            output_mfs = {ouputs.MembershipFunctions};
            
            input1_mf_types = {input_mfs{1}.Type};
            input2_mf_types = {input_mfs{2}.Type};
            output_mf_types = {output_mfs{1}.Type};
            
            
            input1_mf_numericArray = cellfun(@(x) strcmp(x, "trapmf"), input1_mf_types);
            input2_mf_numericArray = cellfun(@(x) strcmp(x, "trapmf"), input2_mf_types);
            output_mf_numericArray = cellfun(@(x) strcmp(x, "trapmf"), output_mf_types);
            
            % MF types array
            mf_types  = [input1_mf_numericArray input2_mf_numericArray output_mf_numericArray];
            
            input1_mf_params = cell2mat({input_mfs{1}.Parameters});
            input2_mf_params = cell2mat({input_mfs{2}.Parameters});
            output_mf_params = cell2mat({output_mfs{1}.Parameters});
            
            % MF parameters array
            mf_params = [input1_mf_params input2_mf_params output_mf_params];
            
            %% Extract the fis fuzzy ranges;
            inputs_range = cell2mat({fis.Inputs.Range}) ;
            output_range = cell2mat({fis.Output.Range});
            
            var_ranges = [inputs_range output_range ];  % fuzzy variable ranges
            
            %% Extract rulebase and Rule weights
            % Get the number of rules
            numRules = length(fis.Rules);
            
            % Initialize an array to hold the rule base in a 2D format
            fis_ruleBase = zeros(numRules, length(fis.Inputs) + length(fis.Outputs));
            
            % Loop through the rules to extract the antecedent and consequent
            for i = 1:numRules
                % Extract antecedents (input conditions)
                antecedents = fis.Rules(i).Antecedent;
                
                % Extract consequents (output actions)
                consequents = fis.Rules(i).Consequent;
                
                % Combine them into a single row for the 2D array
                fis_ruleBase(i, :) = [antecedents consequents];
            end
            rule_base = fis_ruleBase;
            % rule connections
            rule_operator = cell2mat({fis.Rules.Connection});
            rule_connections = rule_operator;
            % Get rule wights as an mat array
            rule_weights = cell2mat({fis.Rules.Weight});
            
        end
        
    end
end