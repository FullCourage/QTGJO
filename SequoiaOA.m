% Sequoia Optimisation Algorithm(SequoiaOA)
% Developed in MATLAB R2022b
function [bestFitness, bestSolution, convergenceCurve] = SequoiaOA(popSize, Max_FEs, lb, ub, dim, fobj)

    % Initial population
    population = repmat(lb, popSize, 1) + rand(popSize, dim) .* (repmat(ub - lb, popSize, 1));

    l=0;
    f=0;
    
    % initializes the fitness array
    fitness = zeros(popSize, 1);
    for i = 1:popSize
        fitness(i) = fobj(population(i, :));
    end
    
    % initializes the best solution
    [bestFitness, bestIndex] = min(fitness);
    bestSolution = population(bestIndex, :);

    f=f+popSize;
    l=l+1;
    convergenceCurve(l)=bestFitness;
    
    % convergence curve
    % convergenceCurve = zeros(1, maxIter);
    
    % Start iteratively updating the solution
    while 1
        % Dynamically adjust parameters
        fireProbability = max(0.3 - 0.15 * (f / Max_FEs), 0.1);
        mutationRate = max(0.2 - 0.1 * (f / Max_FEs), 0.02);
        
        % Sorted by fitness
        [sortedFitness, sortedIndices] = sort(fitness);
        sortedPopulation = population(sortedIndices, :);
        population = sortedPopulation;
        
        eliteSize = 2;
        eliteSolutions = population(sortedIndices(1:eliteSize), :);
        
        % Collective Growth Resource Sharing and Network
        halfPopSize = round(popSize/2);
        topIndividuals = population(1:halfPopSize, :);
        meanTop = mean(topIndividuals, 1);
        population = population + randn(popSize, dim) .* (repmat(meanTop, popSize, 1) - population);
        
        % Adaptation and Resilience
        if rand < fireProbability
            population = population + randn(popSize, dim) * 0.5;
        end
        
        % Reproduction and Diversity
        for i = 1:2:popSize
            if i + 1 <= popSize
                % Crossover
                alpha = rand;
                offspring1 = alpha * population(i, :) + (1 - alpha) * population(i + 1, :);
                offspring2 = alpha * population(i + 1, :) + (1 - alpha) * population(i, :);
                
                % Mutation
                if rand < mutationRate
                    offspring1 = offspring1 + randn(1, dim) * 0.3;
                    offspring2 = offspring2 + randn(1, dim) * 0.3;
                end
                
                % Boundary Treatment
                offspring1 = max(min(offspring1, ub), lb);
                offspring2 = max(min(offspring2, ub), lb);
                
                % Regeneration Population
                population(i, :) = offspring1;
                population(i + 1, :) = offspring2;
            end
        end
        
        % Local Search Mechanism
        localBest = bestSolution + 0.1 * randn(1, dim);
        localBest = max(min(localBest, ub), lb);
        localFitness = fobj(localBest);
        if localFitness < bestFitness
            bestFitness = localFitness;
            bestSolution = localBest;
        end
        f=f+1;
        if f>Max_FEs
            return;
        end
        if mod(f,popSize)==0
            l=l+1;
            convergenceCurve(l)=bestFitness;
        end
        
        % Elite Preservation Mechanism
        population(sortedIndices(end-eliteSize+1:end), :) = eliteSolutions;
        
        % Calculation fitness
        for i = 1:popSize
            fitness(i) = fobj(population(i, :));

            f=f+1;
            if f>Max_FEs
                return;
            end
            if mod(f,popSize)==0
                l=l+1;
                convergenceCurve(l)=bestFitness;
            end
        end
        
        % update the best solution
        [currentBestFitness, currentBestIndex] = min(fitness);
        if currentBestFitness < bestFitness
            bestFitness = currentBestFitness;
            bestSolution = population(currentBestIndex, :);
        end
        
        % Records the convergence curve
        % convergenceCurve(iter) = bestFitness;
       
    end
end