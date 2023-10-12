function out = GAa(problem, params) 

%% Problem 
CostFunction = problem.CostFunction; 
nVar = problem.nVar; 

%% Parameters  
MaxIt = params.MaxIt ;  
Npop = params.nPop ; 
mu = params.mu ; 
betea = params.betea; 
pc = params.pc ; 
nc = (round(pc*Npop)/2)*2; 

SelectionFunc = params.SelectionFunc; 
CrossoverFunc = params.CrossoverFunc; 
% Template for empty individuals 
empty_induvidual.Position = []; 
empty_induvidual.Cost = [] ;  

% Best Genome 
bestGenome.Cost = inf; 

%% Initialization  
population = repmat(empty_induvidual, Npop, 1) ; 

    for i = 1:Npop   
        % Generate random solutions  
        population(i).Position = randi([0,1], 1,nVar); 

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
         popc = repmat(empty_induvidual, nc/2, 2);  

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
             switch CrossoverFunc
                 case "SinglePoint" 
                       [popc(k,1).Position, popc(k,2).Position]= ... 
                            SinglePointCrossover(p1.Position,p2.Position);
                 case "DoublePoint"  
                       [popc(k,1).Position, popc(k,2).Position]= ... 
                            DoublePointCrossover(p1.Position,p2.Position);
                 case "Uniform"  
                       [popc(k,1).Position, popc(k,2).Position]= ... 
                            UniformCrossover(p1.Position,p2.Position);
                 case "Hybrid"  
                       [popc(k,1).Position, popc(k,2).Position]= ... 
                            hybridCrossover(p1.Position,p2.Position);

                 otherwise
                     [popc(k,1).Position, popc(k,2).Position]= ... 
                            hybridCrossover(p1.Position,p2.Position);
             end
         end  % CrossOver
    
        %  Convert popc to Single-Column Matrix
        popc = popc(:); 
    
         % Mutation & Evaluation
         for l=1:nc 
            % Perform mutation
            popc(l).Position = mutation(popc(l).Position, mu); 
            
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

out.pop = population; 
out.bestGenome = bestGenome;  
out.bestcosts = bestcosts; 

%Plot Costs 
% Plot Graph 
figure; 
p = plot(out.bestcosts);
xlabel("Iterations"); 
ylabel("Best Cost"); 
p.Color = 'red'; 
p.LineStyle = '--'; 
p.LineWidth = 1.5; 
grid on ;

end  % end GAa

%% CrossOver functions  

% Single Point cross Over
function [y1, y2] = SinglePointCrossover(x1,x2)    

    nVar = numel(x1) ;
    
    cross_pt = randi([1, nVar-1]);  
    
    y1 = [x1(1:cross_pt) x2(cross_pt+1:end)];  
    y2 = [x2(1:cross_pt) x1(cross_pt+1:end)];  

end 

% Double Point Crossover 
function [y1, y2] = DoublePointCrossover(x1,x2)    

    nVar = numel(x1) ; 
    q = randperm(nVar); 
    cp1 = min(q(1), q(2)); 
    cp2 = max(q(1), q(2)); 

    
    y1 = [x1(1:cp1) x2(cp1+1:cp2) x1(cp2+1:end)]; 
    y2 = [x2(1:cp1) x1(cp1+1:cp2) x2(cp2+1:end)];  

end 

% Uniform CrossOver 
function [y1, y2] = UniformCrossover(x1, x2)
    
    alphea =  randi([0,1], size(x1)); 
    y1 = alphea.*x1 + (1-alphea).*x2;
    y2 = alphea.*x2 + (1-alphea).*x1;

end 

% Random Cross Over Selection  
function [y1, y2] = hybridCrossover(x1, x2) 
    m = randi([1,3]); 
    switch m  
        case 1 
            [y1, y2] = SinglePointCrossover(x1, x2) ;
        case 2 
            [y1, y2] = DoublePointCrossover(x1, x2); 
        case 3 
            [y1, y2] = UniformCrossover(x1, x2); 
        otherwise 
            [y1, y2] = UniformCrossover(x1, x2); 
    end

end

%% Muattion Function 
function y = mutation(x, mu)

    flag = (rand(size(x)) < mu);

    y = x;
    y(flag) = 1 - x(flag);
    
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
