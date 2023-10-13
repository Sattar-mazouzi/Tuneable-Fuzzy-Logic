%% here the problem definition and vraiable decaration are defined 
clear all
fuzzyVarParams= [ 1, 5;  ... 
               5, 5; ... 
               4, 5];  

 mf_tt = randi([0,1], 1,5*3);  
out = genUniformMfs(fuzzyVarParams,mf_tt); 