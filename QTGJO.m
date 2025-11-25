%% QTGJO algorithm
%___________________________________________________________________%
% 1. Spatiotemporal synergy regulation
% 2. Adaptive quantum tunneling
%___________________________________________________________________%

function [ GlobalBestScore, GlobalBest,ConvergenceCurve] = QTGJO(N, Max_FEs, lb, ub, dim, fobj)

    GlobalBest = zeros(1, dim);
    GlobalBestScore = inf;
    SecondBest = zeros(1, dim);
    SecondBestScore = inf;
    
    Positions = initialization(N, dim, ub, lb);
    fitness = inf(N, 1);

    l=0; 
    f=0;

    for i = 1:N
        Positions(i,:) = max(min(Positions(i,:), ub), lb);
        fitness(i) = fobj(Positions(i,:));
        
        if fitness(i) < GlobalBestScore
            SecondBestScore = GlobalBestScore;
            SecondBest = GlobalBest;
            GlobalBestScore = fitness(i);
            GlobalBest = Positions(i,:);
        elseif fitness(i) < SecondBestScore && fitness(i) > GlobalBestScore
            SecondBestScore = fitness(i);
            SecondBest = Positions(i,:);
        end
        f=f+1;
        if f>Max_FEs
            return;
        end
        if mod(f,N)==0
            l=l+1;
            ConvergenceCurve(l)=GlobalBestScore;
        end
    end

    d0=calculateNormalizedDiversity(Positions);

    while 1
        
        E1 = 1.5 * (1 - (f / Max_FEs));
        RL = 0.05 * levy(N, dim, 1.5);

        dt=calculateNormalizedDiversity(Positions);
        if d0 == 0
            dv = 0; 
        else
            dv = dt / d0;
        end

        for i = 1:N
            
            old_pos = Positions(i,:);
            old_fitness = fitness(i);

            mean_pop=mean(Positions);

            if rand<0.5
                for d = 1:dim
                    r1 = rand();
                    E0 = 2*r1 - 1;
                    E = E1 * E0;
                    
                    if dv>0.7
                        M_leader = abs(GlobalBest(d) -  RL(i,d)*Positions(i,d));
                        M(i,d) = GlobalBest(d) - E *M_leader;
                        FM_leader = abs(SecondBest(d) -  RL(i,d)*Positions(i,d));
                        FM(i,d) = SecondBest(d) - E * FM_leader;
                        MEAN_leader = abs(mean_pop(d) -  RL(i,d)*Positions(i,d));
                        MEAN(i,d) = mean_pop(d) - E * MEAN_leader;
                        Positions(i,d) = (M(i,d)+FM(i,d)+MEAN(i,d))/3;
                    elseif dv>0.3
                        M_leader = abs(RL(i,d)*GlobalBest(d) - Positions(i,d));
                        M(i,d) = GlobalBest(d) - E * M_leader;
                        FM_leader = abs(RL(i,d)*SecondBest(d) - Positions(i,d));
                        FM(i,d) = SecondBest(d) - E * FM_leader;
                        
                        Positions(i,d) = (M(i,d)+FM(i,d))/2;
                    else
                        % Rapid convergence
                        M_leader = abs(GlobalBest(d) - Positions(i,d));
                        Positions(i,d) = GlobalBest(d) - rand*E*M_leader + E * (Positions(i,d) - mean_pop(d));
                        
                    end
                end

            else
                for d = 1:dim
                    r1 = rand();
                    E0 = 2*r1 - 1;
                    E = E1 * E0;
                    if abs(E)>1
                        M_leader = abs(GlobalBest(d) -  RL(i,d)*Positions(i,d));
                        M(i,d) = GlobalBest(d) - E *M_leader;
                        FM_leader = abs(SecondBest(d) -  RL(i,d)*Positions(i,d));
                        FM(i,d) = SecondBest(d) - E * FM_leader;
                        MEAN_leader = abs(mean_pop(d) -  RL(i,d)*Positions(i,d));
                        MEAN(i,d) = mean_pop(d) - E * MEAN_leader;
                        Positions(i,d) = (M(i,d)+FM(i,d)+MEAN(i,d))/3;
                    elseif abs(E)>0.5
                        M_leader = abs(RL(i,d)*GlobalBest(d) - Positions(i,d));
                        M(i,d) = GlobalBest(d) - E * M_leader;
                        FM_leader = abs(RL(i,d)*SecondBest(d) - Positions(i,d));
                        FM(i,d) = SecondBest(d) - E * FM_leader;
                        
                        Positions(i,d) = (M(i,d)+FM(i,d))/2;
                    else
                        % Rapid convergence
                        M_leader = abs(GlobalBest(d) - Positions(i,d));
                        Positions(i,d) = GlobalBest(d) - rand*E*M_leader + E * (Positions(i,d) - mean_pop(d));
                        
                    end
    
                end
            end
            
            Positions(i,:) = max(min(Positions(i,:), ub), lb);
            new_fitness = fobj(Positions(i,:));

            f=f+1;
            if f>Max_FEs
                return;
            end
            if mod(f,N)==0
                l=l+1;
                ConvergenceCurve(l)=GlobalBestScore;
            end

            if new_fitness < old_fitness
                fitness(i) = new_fitness;
                if new_fitness < GlobalBestScore
                    SecondBestScore = GlobalBestScore;
                    SecondBest = GlobalBest;
                    GlobalBestScore = new_fitness;
                    GlobalBest = Positions(i,:);
                end
            else
                Positions(i,:) = old_pos; 
                fitness(i) = old_fitness;
            end
            old_pos = Positions(i,:);
            old_fitness = fitness(i);

            if GlobalBestScore == 0
                fitness_diff = abs(old_fitness - GlobalBestScore); 
            else
                fitness_diff = (old_fitness - GlobalBestScore) / abs(GlobalBestScore);
            end
            barrier_height = tanh(10* fitness_diff);

            
            % Quantum tunneling probability
            tunneling_prob =exp(-barrier_height / (0.1 + 0.9*(1 - f/Max_FEs))); 
            

            if rand() < tunneling_prob
                % Tunneling successful

                rand_dim = randi(dim);
                scale_ratio = (Positions(i, rand_dim) - lb(rand_dim)) / (ub(rand_dim) - lb(rand_dim));
                new_pos = lb + scale_ratio * (ub - lb);

                direction = randn(1,dim);
                direction = direction / norm(direction);

                Positions(i, :) = new_pos + direction .* norm(Positions(i, :) - new_pos)*rand^2;

            else
                % Tunneling failure
                R = randperm(N, 2);
                while any(R == i)
                    R = randperm(N, 2);
                end

                triangle_center = (Positions(i,:) + Positions(R(1),:)  + Positions(R(2),:) ) / 3;
                direction = randn(1,dim);
                direction = direction / norm(direction);
                Positions(i,:) = triangle_center +  direction * (1 - f / Max_FEs)^0.5.*norm(triangle_center-Positions(i,:));
            end

            Positions(i,:) = max(min(Positions(i,:), ub), lb);

            new_fitness = fobj(Positions(i,:));

            if new_fitness < old_fitness
                
                fitness(i) = new_fitness;

                if new_fitness < GlobalBestScore
                    SecondBestScore = GlobalBestScore;
                    SecondBest = GlobalBest;
                    GlobalBestScore = new_fitness;
                    GlobalBest = Positions(i,:);
                end
            else
                Positions(i,:) = old_pos; 
                fitness(i) = old_fitness;
            end

            f=f+1;
            if f>Max_FEs
                return;
            end
            if mod(f,N)==0
                l=l+1;
                ConvergenceCurve(l)=GlobalBestScore;
            end
            
        end
       
    end
end

function Positions = initialization(N, dim, ub, lb)
    if numel(ub) == 1
        Positions = rand(N, dim) .* (ub - lb) + lb;
    else
        Positions = zeros(N, dim);
        for d = 1:dim
            Positions(:, d) = rand(N, 1) .* (ub(d) - lb(d)) + lb(d);
        end
    end
end

%% Population diversity
function dt = calculateNormalizedDiversity(Positions)
    
    current_center = mean(Positions, 1);  
    
    current_distances = sqrt(sum((Positions - current_center).^2, 2));
    dt = mean(current_distances);
end