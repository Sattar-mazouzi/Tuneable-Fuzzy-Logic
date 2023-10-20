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
