# Tuneable Fuzzy Logic

_This project is functional but still under development to improve functionality and add more features._

## Overview
![Alt text](Tunable_flc.png)
The Tunable Mamdani Fuzzy Logic Inference System (M-FIS)  is a customizable and optimizable fuzzy logic system. This system allows users to fine-tune their FIS through various optimization algorithms to meet specific cost function requirements. Users can create up to two inputs and one output, with the ability to customize each fuzzy variable using triangular or trapezoidal membership functions.

## Features

- **Flexible FIS Creation**: Easily create a two-input, one-output FIS.
- **Optimization Algorithms**: Choose from different optimization algorithms to fine-tune your FIS.
- **Fuzzy Variable Range Optimization**: Optimize the ranges of fuzzy variables.
- **Membership Function Types**: Customize membership functions by choosing between triangular and trapezoidal types.
- **Membership Function Parameters**: Tune the parameters and shapes of membership functions.
- **Rule-Base Tuning**: Optimize the rule-base for better performance.
- **Rule-Weights Optimization**: Fine-tune the weights of the rules.
- **Comperhensive Optimization**: Optimized the parameters at once (limited in this version it's possible to optimzed rule-base rule weights and memebership parameters).

## Usage Manual

### 1. Define the Objective Function

Specify the objective function that needs to be minimized. This function will guide the optimization process.
the objective function should take the following form:

```matlab
function out = fitness(fis,X)  
    fis.evaluate(X); 

    cost  =  xxx

out = cost ;
end
```

Where X: is the decision variable, Cost: is the objective value that can be calculated after evaluating the fis to new parameters.

### 2. Create the FIS

Initialize the FIS with the desired parameters:

```matlab
% input 1 name, fuzzy range and  the number of membership function 
fis_in1 = struct(...); % Define input 1 parameters
fis_in1.name = "in1_name" ;      
fis_in1.range = [lower upper];         
fis_in1.MfNumber = 3; 
fis_in1.MfTypes = [mf1_type mf2_type mf3_type]; %0 for Triangular ,1 for trapezoidal  

% input 2 name, fuzzy range and  the number of membership function
fis_in2 = struct(...); % Define input 2 parameters
fis_in2.name = "in2_name" ;      
fis_in2.range = [lower upper];         
fis_in2.MfNumber = 3; 
fis_in2.MfTypes = [mf1_type mf2_type mf3_type];

% input 1 name, fuzzy range and  the number of membership function
fis_out = struct(...); % Define output parameters
fis_out.name = "out_name" ;      
fis_out.range = [lower upper];         
fis_out.MfNumber = 3; 
fis_out.MfTypes = [mf1_type mf2_type mf3_type];

fis_rules = [...]; % Define the rule base

tuneableFis = TuneableFis(fis_in1, fis_in2, fis_out, fis_rules); 

% generate a uniform distributed mfs fis   
fis = tuneableFis.Fis 
```

_**Note:**_ After each optimization phase, the Fis will be updated to the new  optimized Fis.

### 3. Define optimization parameters

Define a struct that holds the optimization parameters:

```matlab
problem = struct(...); % Define the optimization problem parameters 
problem.CostFunction = @(x) fitness(tuneableFis,x); %the fitness (objective) function 
problem.MaxIt= 5;   % number of iteration 
problem.nPop = 10;  % number of population
```

### 4. Optimize Membership Function Types

Optimize the types of membership functions:

```matlab
output = tuneableFis.TuneMfTypes(problem);
```

### 5. Optimize Membership Function Parameters

Tune the parameters and shapes of the membership functions:

```matlab
output = tuneableFis.TuneMfParameters(problem);
```

### 6. Optimize Rule-Base

Optimize the rule-base for better performance:

```matlab
output = tuneableFis.TuneRules(problem);
```

In tune Rules, the TuneRules function can take other tow parameters tune_operator, rule_limits.

- **tune_operator**: which is a boolean type that can be either **true** or **false**  indecating tune/not tune the AND/OR connections in the rules.

- **rule_limits**: is 1-D array specifying the limites of the output subsets for each rule.

### 7. Optimize Rule Weights

Fine-tune the weights of the rules:

```matlab
output = tuneableFis.TuneRuleWeigths(problem);
```

### 8. Optimize Fuzzy Variable Ranges

Optimize the ranges of fuzzy variables:

```matlab
output = tuneableFis.TuneFuzzyVarsRange(problem);
```

The TuneFuzzyVarsRange method take two other parameter:

- **tune_range** : which has three values {"upper", "lower", "both"} specifying, "upper" to tune the upper bound, "lower" to tune the lower bound and the value "both" for tuning both sides.  
- **range_offset**: a number between 0.1 and 2 specifying how much to tune the fuzzy range sides.

### 9. Comperhensive Optimization

Optimize the MF parameters, rule-base and rule weights all at once.

```matlab
output = Optimizefis(problem);
```

## Conclusion

This project represents a significant advancement in the field of fuzzy logic, providing a powerful and flexible tool for optimizing complex systems. The Tuneable FIS system is a valuable asset for researchers, engineers, and anyone looking to optimize their processes and achieve better outcomes.

## Future Work

* **Enhanced Functionality:** Continue to improve and add more features.

- **User Interface:** Develop a user-friendly interface for easier interaction.

- **Documentation:** Provide comprehensive documentation and examples.

Thank you for using the Tuneable Fuzzy Logic Inference System. We hope it proves to be a valuable tool in your optimization endeavors.
