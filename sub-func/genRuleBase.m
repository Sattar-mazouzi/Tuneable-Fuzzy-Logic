function out = genRuleBase(fisVarParameters,rules_data)

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