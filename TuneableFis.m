classdef  TuneableFis
    properties
     input1 ;  
     input2; 
     output;
     Rule_Base; 
     Fis;  
    end 
    methods
        function obj = TuneableFis(fis_in1, fis_in2, fis_out, fis_rules)
            if nargin > 0 
                obj.input1 = fis_in1; 
                obj.input2 = fis_in2; 
                obj.output = fis_out; 
                obj.Rule_Base = fis_rules;
            else
                disp("you must specify the Fis parameters"); 
            end
        end
        function fis = create_fis(obj)
            in1 = obj.input1;
            in2 =  obj.input2; 
            out = obj.output;
            rules = obj.Rule_Base; 
            mftypes = {[in1.MfTypes], [in2.MfTypes], [out.MfTypes]}; 

            fis = tunebale_flc(obj.input1, obj.input2, obj.output,mftypes, rules,0) ;
        
        end
        function tuned_fis = tune_fis(obj)
            in1 = obj.input1;
            in2 =  obj.input2; 
            out = obj.output;
            % Fuzzy variable paramters 
            fuzzyVarParams = [in1.range, in1.MfNumber;   
            in2.range, in2.MfNumber;
            out.range, out.MfNumber] ;  
            

            tuned_fis = 0; 
        end

    end 
end