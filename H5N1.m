%_________________________________________________________________________________
%  H5N1 algorithm source codes version 1.0
%
%  Developed in MATLAB R2022a
%
%  Author and programmer: Le Xuan Thang
%
%         e-Mail: lexuanthang.official@gmail.com
%                 lexuanthang.official@outlook.com
%         website: https://lexuanthang.vn
%
%   Main paper: Le, T.X., Bui, T.T. and Tran, H.N. (2025), "The H5N1 algorithm: a viral-inspired optimization for solving real-world engineering problems", Engineering Computations 
%   DOI:https://doi.org/10.1108/EC-05-2024-0472
%____________________________________________________________________________________
function [GlobalBest_Cost, GlobalBest_Position, BestCosts] = H5N1(nvirus, Max_FEs, lb, ub, dim, fobj)

if size(lb,2) == 1
    lb = repmat(lb,1,dim);
    ub = repmat(ub,1,dim);
end

% setting up parameter for H5N1 algorithms
P1 = 0.8; % Probability of attaching to a bird, poultry
P2 = 0.85; % Probability of survival after infection

w = 1;
wmin = 0.001;
wmax = 1;
p = 0;
m_virus = zeros(nvirus,dim);                     % mutate for virus virus population
m_virusBest = zeros(nvirus,dim);                 % mutate for virusBest virus population
ro = (0:1:nvirus-1);                             % rotating index array
R = rand(nvirus,dim);                              % matrix of normalization random number for all
v = rand(nvirus,1);                                % vector of normalization random number for mutate

%% Initialize the virus population randomly
virusCost = zeros(nvirus,1);                       % cost of virus of virus
virusPosition  = zeros(nvirus,dim);                % position of virus of virus

virusBestPosition = rand(nvirus, dim).*(ub-lb) + lb; % initilization first virus population
virusBestCost = zeros(nvirus, 1);


l=0; 
f=0;
% Evaluate the objective function for each individual
for i = 1:nvirus
    virusBestCost(i) = fobj(virusBestPosition(i, :));
end

% Find the best individual
[GlobalBest_Cost, bestidx] = min(virusBestCost);
GlobalBest_Position = virusBestPosition(bestidx, :);

f=f+nvirus;
if mod(f,nvirus)==0
    l=l+1;
    BestCosts(l) = GlobalBest_Cost;
end

%% value store for visualization
% BestCosts = zeros(Max_FEs,1);                    % save best solution of every iteration

%% Main loop
while 1
    c = exp(-4*f/Max_FEs)*rand(); % Eq. 3-9
    p = 1 ./ (1 + exp(-10*(f/Max_FEs-0.5))) + rand()*(1-p); % Eq. 3-10 Multation probability
    wd = wmin + (wmax-wmin)*exp(-(f/Max_FEs)); % Eq. 3-8 (part 2)

    % Determine the probability of attaching to a bird
    P_attack = rand(nvirus, 1);
    % Determine the probability of survival after infection
    P_adapt = rand(nvirus, 1);
    
    sigma = randperm(3);                 % index pointer array
    a1 = randperm(nvirus);               % Eq. 3-1
    rt1 = rem(ro+sigma(1),nvirus);       % Eq. 3-2
    a2 = a1(rt1+1);                      % Eq. 3-3
    rt2 = rem(ro+sigma(2),nvirus);       % Eq. 3-4
    a3 = a2(rt2+1);                      % Eq. 3-5

    pmp1 = virusBestPosition(a1,:);    % permutate virus population 1
    pmp2 = virusBestPosition(a2,:);    % permutate virus population 2
    pmp3 = virusBestPosition(a3,:);    % permutate virus population 3

    v  = v(a3,:);                 % permutate coefficient
    sm = rand(nvirus,dim);        % stochastic for mutate

    % Update the virus population using the infected fitness
    for i = 1:nvirus
        m_virus(i,:) = sm (i,:) < v(i,1);                    % all random numbers < vm are 1, 0 otherwise
        m_virusBest(i,:) = m_virus (i,:) < v(i,1);           % inverse mutate virus to virusBest

        % Move the individual according to the probability of survival
        % Check if the virus attaches to poultry or human
        if P_attack(i) < P1 % attach to poultry

            % Calculate the probability of survival of the poultry
            if P_adapt(i) < P2 % Good enviroment / the poultry will survives
                % Get direction
                direction = GlobalBest_Position - virusBestPosition(i,:);
                virusPosition(i,:) = pmp3(i,:)+ R(i,:).*direction; % Eq. 3-6 (part 1)
            
            else % Bad enviroment / the poultry will not survive
                xold = GlobalBest_Position ;                                                    % Eq. 3-7 (part 1)
                xnew = virusPosition(i,:) + R(i,:).*(GlobalBest_Position-virusPosition(i,:));   % Eq. 3-7 (part 2)
                virusPosition(i,:) = 1/2*c*w*(xnew+xold)*rand();                                % Eq. 3-7 (part 3)
                
                virusPosition(i, :) = max(min(virusPosition(i, :), ub), lb);
                virusCost(i) = fobj(virusPosition(i,:));

                ps = rand();
                % Update the global best position and cost Eq. 3-12
                if ps < p && virusCost(i) < GlobalBest_Cost
                    GlobalBest_Position = virusPosition(i,:);
                    GlobalBest_Cost = virusCost(i);
                end
                f=f+1;
                if f>Max_FEs
                    return;
                end
                if mod(f,nvirus)==0
                    l=l+1;
                    BestCosts(l) = GlobalBest_Cost;
                end
            end

        else % attack to human
            if P_adapt(i) < P2 % Good enviroment /the human will survives
                % Get direction
                direction = GlobalBest_Position - virusBestPosition(i,:);
                virusPosition(i,:) = virusPosition(i,:)  + (R(i,:).*(pmp1(i,:)-pmp2(i,:)) + R(i,:).*direction)/2; % Eq. 3-6 (part 2)

            else % Bad enviroment / the human will not survive
                xold = GlobalBest_Position ;                                                        % Eq. 3-7 (part 1)
                xnew = virusPosition(i,:) + R(i,:).*(GlobalBest_Position-virusBestPosition(i,:));   % Eq. 3-7 (part 2)
                virusPosition(i,:) = 1/2*c*w*(xnew+xold)*rand();                                    % Eq. 3-7 (part 3)

                virusPosition(i, :) = max(min(virusPosition(i, :), ub), lb);
                virusCost(i) = fobj(virusPosition(i,:));

                ps = rand();
                % Update the global best position and cost Eq. 3-12
                if ps < p && virusCost(i) < GlobalBest_Cost
                    GlobalBest_Position = virusPosition(i,:);
                    GlobalBest_Cost = virusCost(i);
                end
                f=f+1;
                if f>Max_FEs
                    return;
                end
                if mod(f,nvirus)==0
                    l=l+1;
                    BestCosts(l) = GlobalBest_Cost;
                end

            end
        end
    end

    % Mutate
    virusPosition = virusPosition.* m_virus + virusBestPosition.* m_virusBest; % Eq. 3-11

    %% CheckBoundary / Clip the position to the search space
    changeRows =  virusPosition > ub;
    virusPosition(changeRows) =  pmp1(changeRows);
    changeRows =  virusPosition < lb;
    virusPosition(changeRows) = pmp2(changeRows);

    % Evaluate new position
    for i = 1:nvirus
        virusCost(i,:) = fobj(virusPosition(i,:));   % check cost of competitor

        if (virusCost(i) <= virusBestCost(i))
            virusBestPosition(i,:) = virusPosition(i,:);
            virusBestCost(i)   = virusCost(i);
            if (virusBestCost(i) <= GlobalBest_Cost)
                GlobalBest_Cost = virusBestCost(i);
                GlobalBest_Position = virusPosition(i,:);
            end
        end

        f=f+1;
        if f>Max_FEs
            return;
        end
        if mod(f,nvirus)==0
            l=l+1;
            BestCosts(l) = GlobalBest_Cost;
        end
    end

    % Update for virusBest and Globalbest Eq. 3-12
%     for i = 1:nvirus
%         if (virusCost(i) <= virusBestCost(i))
%             virusBestPosition(i,:) = virusPosition(i,:);
%             virusBestCost(i)   = virusCost(i);
%             if (virusBestCost(i) <= GlobalBest_Cost)
%                 GlobalBest_Cost = virusBestCost(i);
%                 GlobalBest_Position = virusPosition(i,:);
%             end
%         end
%     end

    w = w * wd; % Eq. 3-8 (part 1)
%     BestCosts(iter) = GlobalBest_Cost;
%     if (mod(iter, 1) == 0)
%         fprintf('The best solutions at iteration %4.0f is: %3.12f \n',iter, GlobalBest_Cost)
%     end

end

