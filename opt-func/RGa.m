function out = RGA(problem, params) 
 
%% Problem 
CostFunction = problem.CostFunction; 
nVar = problem.nVar;
VarMin = problem.VarMin;
VarMax = problem.VarMax; 
VarSize = [1 nVar] ; 

%% Parameters  
MaxIt = params.MaxIt ;  
Npop = params.nPop ; 
mu = params.mu ; 
sigma = params.sigma; 
betea = params.betea; 
gamma = params.gamma ;
pc = params.pc ; 
nc = (round(pc*Npop/2))*2; 


SelectionFunc = params.SelectionFunc; 
% Template for empty individuals 
empty_induvidual.Position = []; 
empty_induvidual.Cost = [] ;  

% Best Genome 
bestGenome.Cost = inf; 

%% Initialization  

population = repmat(empty_induvidual, Npop, 1) ; 

    for i = 1:Npop  
       
        % Generate random solutions  
        population(i).Position = unifrnd(VarMin, VarMax, VarSize); 
        % Evaluate Solution
        population(i).Cost = CostFunction(population(i).Position);

        % update the bestsolution  
        if population(i).Cost < bestGenome.Cost
            bestGenome = population(i); 
        end
    end 

% Iteration Bets cost Histtory  
bestcosts = nan(MaxIt,1); 

    %% Main Loop 
    for it = 1:MaxIt
    
         %Initialize offsprings Population 
         popc = repmat(empty_induvidual,nc/2, 2);  

         % Selections Probabilities  
         probs = sel_probability(population, betea);
    
         %CrossOver 
         for k =  1:nc/2 
             
             %Selcet parents  % and generate offsprings  
             switch SelectionFunc
                 case "Random" 
                     idx = random_select(Npop); 
                 case "RuletteWheel"
                     idx = RuletteWheelSelection(probs); 
                 otherwise
                     idx = hybridSelection(Npop, probs); 
             end
             p1 = population(idx); 
             p2 = population(idx);
             
             % porform CrossOver 
              [popc(k,1).Position, popc(k,2).Position]= ... 
                       UniformCrossover(p1.Position,p2.Position, gamma);

         end  % CrossOver
    
        %  Convert popc to Single-Column Matrix
        popc = popc(:); 
    
         % Mutation & Evaluation
         for l=1:nc 

            % Perform mutation
            popc(l).Position = mutation(popc(l).Position, mu, sigma);  

            % Check for bounds  
            popc(l).Position =  max(popc(l).Position, VarMin); 
            popc(l).Position =  min(popc(l).Position, VarMax); 

            
            % Evaluation
            popc(l).Cost = CostFunction(popc(l).Position);
     
            % update the bestsolution  
            if popc(l).Cost < bestGenome.Cost
                bestGenome = popc(l); 
            end
         end % Mutation 
    
 
    
        % Merge Population  
        population = [population; popc]; 
    
        % Sort population 
        population = sort_pop(population); 
          
        % Remove Extra Individuals 
        population = population(1:Npop); 

        % Update Best costof iterations 
        bestcosts(it) = bestGenome.Cost; 
        
        % Display Infos 
        disp(['Iteration ' num2str(it) ': Best Cost = ' ...
            num2str(bestcosts(it))]); 
      
    
    
    end % Main Loop

% Outputs 
out.pop = population; 
out.bestGenome = bestGenome;  
out.bestcosts = bestcosts; 

% Plot Graph 
figure; 
p = semilogy(out.bestcosts);
xlabel("Iterations"); 
ylabel("Best Cost"); 
p.Color = 'red'; 
p.LineStyle = '--'; 
p.LineWidth = 1.5; 
grid on ;


end  % end RGA

%% CrossOver functions  

% Uniform CrossOver 
function [y1, y2] = UniformCrossover(x1, x2, gamma)
    
%   alphea = rand(size(x1));  
    alphea = unifrnd(-gamma, 1+gamma, size(x1)); 
    y1 = alphea.*x1 + (1-alphea).*x2;
    y2 = alphea.*x2 + (1-alphea).*x1;

end 



%% Muattion Function 
function y = mutation(x, mu, sigma)

    flag = (rand(size(x)) < mu);
    y = x;
    r = randn(size(x));   
    y(flag) = x(flag) + sigma*r(flag) ;
    
end

%% sort 
function out = sort_pop(pop)
    [~, sor] = sort([pop.Cost]); 
    out = pop(sor); 
end 

%% Selections  

% Random selection 
function idx = random_select(npop)  

        randn = randperm(npop); 
        r = randi([1,npop]);
        idx = randn(r); 

end

% Rullete wheel selection 
function idx = RuletteWheelSelection(probs)
    r = rand*sum(probs); 
    c = cumsum(probs); 
    idx = find(r<=c, 1, 'first'); 
end

% Hybrid Selection 
function idx = hybridSelection(npop,probs) 
    rr = randi([1,2]);  
    switch rr 
        case 1 
            idx = RuletteWheelSelection(probs) ;
        case 2 
            idx = random_select(npop) ; 
        otherwise 
            idx = RuletteWheelSelection(probs); 
    
    end 
end

%% Selection Probabilities 

function probs = sel_probability(pop, betea)
    c = [pop.Cost] ; 
    avgc = mean(c);
    if avgc ~= 0
        c= c/avgc; 
    end
    probs = exp(-betea*c); 
end 