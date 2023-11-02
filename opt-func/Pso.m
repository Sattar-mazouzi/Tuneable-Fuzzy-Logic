function [out,logg, hist] = Pso(problem, params ) 


    %% Problem definistion 
    
    CostFunction = problem.CostFunction;  % Cost Function 
    nVar = problem.nVar ;           % Number Of Unkown Variable 
    VarSize = [1 nVar]  ;           % Matrix Size of Decision Variables  
    VarMin = problem.VarMin ;       % Lower Boun of the Variables 
    VarMax = problem.VarMax  ;      % Upper Bound of Decision Variable 

    %low_bnd = problem.lower_bd ; 
    %up_bnd = problem.upper_bd; 
    
    %% Parameters of PSO 
    stopping = params.stop; 
    MaxIt = params.MaxIt;  % Maxinum Number Of Iteration 
    nPop =  params.nPop  ;    % Population Size  
    W  =  params.W ;        % Inertia Coefficient 
    C1 = params.C1 ;         % Personal Coefficient 
    C2 = params.C2  ;       % Social Coefficient 
    Wdamp =  params.Wdamp ; % damping coefficient of the inertia wheight 
    
    MaxVel = 0.3*(VarMax - VarMin); 
    MinVel = - MaxVel; 
    tolerance = problem.tolerance; 
    max_conv = problem.max_conv; 

    %% Initialization 
    it = 1;  % initialize the iteration count. 
    convergence_count = 0; %initialize the convergence counter
    empty_particle.Position = [] ; 
    empty_particle.Velocity = [] ;
    empty_particle.Cost = [];
    empty_particle.Best.Position = [] ;
    empty_particle.Best.Cost = [];  
    
    % Creat array of empty particls (Population)
    Particle = repmat(empty_particle, nPop, 1); 
    
    % Initialize the Global best 
    GlobalBest.Cost = inf; 
    
    
    for i = 1:nPop 
          % Generate Random Solution
         Particle(i).Position = unifrnd(VarMin, VarMax, VarSize);
         %Particle(i).Position = randi(10*[VarMin VarMax], VarSize)/10;

         % Set the constraint of MFs
         %pose = Particle(i).Position ; 

         %low_pose = [-1 -1 pose(1) pose(3) pose(2) pose(5) pose(4) pose(7)];  
         %up_pose  = [pose(3) pose(5) pose(4) pose(7) pose(6) 1 pose(8) 1];  
         
         %% aapply bounds 
         Particle(i).Position = max(Particle(i).Position, VarMin); 
         Particle(i).Position = min(Particle(i).Position, VarMax);
    
         % Initialize Velocity 
         Particle(i).Velocity = zeros(VarSize); 
        
         % Evaluation 
         Particle(i).Cost = CostFunction(Particle(i).Position);  
    
         % Update the Poersonal Best 
         Particle(i).Best.Position = Particle(i).Position; 
         Particle(i).Best.Cost = Particle(i).Cost;
    
         % Update the Global Best 
         if Particle(i).Best.Cost < GlobalBest.Cost 
             GlobalBest = Particle(i).Best;  
         end
    end 
    
    BestCosts = zeros(MaxIt,1) ;
    
    %% Main Loop Of PSO 
   % for it=1:MaxIt 
              while (convergence_count < max_conv  && it <= MaxIt)  
    
                for i = 1:nPop 
        
                            % Update Velocity 
                        Particle(i).Velocity = W*Particle(i).Velocity ... 
                                           + C1*rand(VarSize).*(Particle(i).Best.Position - Particle(i).Position) ...
                                           + C2*rand(VarSize).*(GlobalBest.Position -Particle(i).Position);  
    
                        % Limit Velociy  
                                
                          Particle(i).Velocity = max(Particle(i).Velocity, MinVel); 
                        Particle(i).Velocity = min(Particle(i).Velocity, MaxVel);
            
                            % Update Position 
                        Particle(i).Position = Particle(i).Position + Particle(i).Velocity; 
    
                           % Set the constraint of MFs
                      
                        % applly lower bound and upper bound 
    
                        Particle(i).Position = max(Particle(i).Position, VarMin); 
                        Particle(i).Position = min(Particle(i).Position, VarMax);
    
                        % Evaluation 
                        Particle(i).Cost = CostFunction(Particle(i).Position); 
                        
    
    
                        if Particle(i).Cost < Particle(i).Best.Cost 
        
                            % update personal best 
                            Particle(i).Best.Cost  = Particle(i).Cost ; 
                            Particle(i).Best.Position = Particle(i).Position ; 
                                
                            % Update the global best 
                             if Particle(i).Cost < GlobalBest.Cost 
                                
                                 GlobalBest = Particle(i).Best; 
                                 
                             end 
                    
                    end 
                
                disp(['it: ' num2str(it) ' particle: ' num2str(i)  ' BestCost = ' num2str(Particle(i).Best.Cost)]);
                end
   
            best_hist(it) = GlobalBest;
            BestCosts(it) = GlobalBest.Cost; 
          % stoppin critirea ;  
         if stopping == true
         
               % check if the best particle's position has not changed for a certain number of consecutive iterations
            if it>1 
                dist = norm(GlobalBest.Position - best_hist(it-1).Position,2);
                disp(dist); 
                if dist < tolerance 
                    convergence_count = convergence_count +1; 
                else 
                    convergence_count = 0; 
                end 
            end 
         end 
           
            % damping inertia coeff 
            W = W*Wdamp;
            disp(['Iteration Results: ' num2str(it) ]) 
            disp(GlobalBest)
         it = it +1;
       end
    %end




out = GlobalBest ;
logg = BestCosts;  
hist = best_hist; 
%% Results 
%Plot the best costs 
figure; 
semilogy(logg,'LineWidth',2,'LineStyle', '-.','Color', 'r'); 
xlabel('Iteration'); ylabel('Best Cost');  
grid on; 

end 

