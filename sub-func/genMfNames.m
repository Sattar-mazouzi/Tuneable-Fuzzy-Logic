function out = genMfNames(mftypes) 
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