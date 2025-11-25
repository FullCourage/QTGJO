%% Q-learning
function [Male_jackal_fit, Male_jackal, Convergence] = QLGJO(N, Max_FEs, lb, ub,  dim, fobj)
%% Initialize parameters
    Male_jackal = zeros(1, dim);
    Male_jackal_fit = inf;

    Female_jackal = zeros(1, dim);
    Female_jackal_fit = inf;

    % initialize the params of Q-Learning
    action_num = 3;
    Reward_table = zeros(action_num, action_num, N);
    Q_table = zeros(action_num, action_num, N);
    cur_state = randi(action_num);

    gamma = 0.5;
    lambda_initial = 0.9;
    lambda_final = 0.1;

    x=initalization(N,dim,ub,lb);

    l=0; 
    f=0;
    for i=1:N
        Flag4ub=x(i,:)>ub;
        Flag4lb=x(i,:)<lb;
        x(i,:)=(x(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        fitness(i)=fobj(x(i,:));

        if fitness(i)<Male_jackal_fit
            Male_jackal_fit=fitness(i);
            Male_jackal=x(i,:);
        end

        if fitness(i)>Male_jackal_fit && fitness(i) < Female_jackal_fit
            Female_jackal_fit=fitness(i);
            Female_Jackal=x(i,:);
        end

        f=f+1;
        if f>Max_FEs
            return;
        end
        if mod(f,N)==0
            l=l+1;
            Convergence(l)=Male_jackal_fit;
        end
    end

    %% Iteration
    while 1
        rl = 0.05 * levy(N, dim, 1.5);
        E1 = 1.5 * (1 - (f / Max_FEs));

        

        for i = 1:N
            for j = 1:dim
                % r1 is a random number in [0, 1]
                r1 = rand();
                E0 = 2 * r1 - 1;
                % Evading energy
                E = E1 * E0;

                if (Q_table(cur_state, 1, i) >= Q_table(cur_state, 2, i) && Q_table(cur_state, 1, i) >= Q_table(cur_state, 3, i))
                    action = 1;
                    % Exploration
                    D_male_jackal = abs((Male_jackal(j) - rl(i, j) * x(i, j)));
                    Male_Positions(i, j) = Male_jackal(j) - E * D_male_jackal;
                    D_female_jackal = abs((Female_jackal(j) - rl(i,j) * x(i, j)));
                    Female_Positions(i, j) = Female_jackal(j) - E * D_female_jackal;
                elseif  (Q_table(cur_state, 2, i) >= Q_table(cur_state, 1, i) && Q_table(cur_state, 2, i) >= Q_table(cur_state, 3, i))
                    action = 2;
                    if rand > 0.5
                        D_male_jackal = abs((Male_jackal(j) - rl(i, j) * x(i, j)));
                        Male_Positions(i, j) = Male_jackal(j) - E * D_male_jackal;
                        D_female_jackal = abs((Female_jackal(j) - rl(i,j) * x(i, j)));
                        Female_Positions(i, j) = Female_jackal(j) - E * D_female_jackal;
                    else
                        D_male_jackal = abs((rl(i, j) * Male_jackal(j) - x(i, j)));
                        Male_Positions(i, j) = Male_jackal(j) - E * D_male_jackal;
                        D_female_jackal = abs((rl(i, j) * Female_jackal(j) - x(i, j)));
                        Female_Positions(i, j) = Female_jackal(j) - E * D_female_jackal;
                    end
                else
                    action = 3;
                    % Exploitation
                    D_male_jackal = abs((rl(i, j) * Male_jackal(j) - x(i, j)));
                    Male_Positions(i, j) = Male_jackal(j) - E * D_male_jackal;
                    D_female_jackal = abs((rl(i, j) * Female_jackal(j) - x(i, j)));
                    Female_Positions(i, j) = Female_jackal(j) - E * D_female_jackal;
                end
                Position(j) = Male_Positions(i, j) + Female_Positions(i, j) / 2;
            end
            Flag4lb = Position < lb;
            Flag4ub = Position > ub;
            Position = Position .* (~(Flag4ub + Flag4lb)) + ub .* Flag4ub + lb .* Flag4lb;
            new_fitness = fobj(Position);

            f=f+1;
            if f>Max_FEs
                return;
            end
            if mod(f,N)==0
                l=l+1;
                Convergence(l)=Male_jackal_fit;
            end

            if action == 1
                % exoplore mutation
                mutation = x(randi(N), :) + 0.85 .* (x(randi(N), :) - x(randi(N), :)) + 0.85 .* (x(randi(N), :) - x(randi(N), :));
            elseif action == 2
                % hybird mutation
                mutation = x(randi(N), :) + 0.85 .* (x(randi(N), :) - x(randi(N), :));
            else
                % exploitation mutation
                mutation = x(i, :) + 0.85 .* (x(randi(N), :) - x(randi(N), :));
            end

            Flag4lb = mutation < lb;
            Flag4ub = mutation > ub;
            mutation = mutation .* (~(Flag4ub + Flag4lb)) + ub .* Flag4ub + lb .* Flag4lb;
            mutation_fit = fobj(mutation);
            if mutation_fit < new_fitness
                new_fitness = mutation_fit;
                Position = mutation;
            end

            if new_fitness < fitness(i)
                x(i, :) = Position;
                fitness(i) = new_fitness;
                Reward_table(cur_state, action, i) = 1;
            else
                Reward_table(cur_state, action, i) = -1;
            end

            if fitness(i)<Male_jackal_fit
                Male_jackal_fit=fitness(i);
                Male_jackal=x(i,:);
            end
    
            if fitness(i)>Male_jackal_fit && fitness(i) < Female_jackal_fit
                Female_jackal_fit=fitness(i);
                Female_Jackal=x(i,:);
            end
    
            f=f+1;
            if f>Max_FEs
                return;
            end
            if mod(f,N)==0
                l=l+1;
                Convergence(l)=Male_jackal_fit;
            end

            % update the Q_table
            r =  Reward_table(cur_state, action, i);
            maxQ = max(Q_table(action, :, i));
            lambda = (lambda_initial + lambda_final) / 2 - (lambda_initial - lambda_final) / 2 * cos(pi * (1 - f / Max_FEs));
            Q_table(cur_state, action, i) = Q_table(cur_state, action, i) + lambda * (r + gamma * maxQ - Q_table(cur_state, action, i));
            cur_state = action;
        end

    end
end


function X=initalization(SearchAgents_no,dim,ub,lb)
    Boundary_no=size(ub,2); 

    if Boundary_no==1 
        X=rand(SearchAgents_no,dim).*(ub-lb)+lb; 
    end

    if Boundary_no>1 
        for i=1:dim
            ub_i=ub(i);
            lb_i=lb(i);
            X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
        end
    end

end
