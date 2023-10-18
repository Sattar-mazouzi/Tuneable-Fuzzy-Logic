function [out1, out2 ]= genUniformMfs(fuzzyVarParams, mf_types)
mfs = {[],[],[]}; 
var1  = fuzzyVarParams(1,:); 
var2 = fuzzyVarParams(2,:);
var3 = fuzzyVarParams(3,:);
mf_types_cell = {[mf_types(1:var1(3))],[mf_types(var1(3)+1:var1(3)+var2(3))], ...
    [mf_types(var1(3)+var2(3)+1:var1(3)+var2(3)+var3(3))]} ;
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
out2 = mf_types_cell; 
end 