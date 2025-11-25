% DE
function [S,P,Convergence_curve]  = DE(nPop,Max_FEs,VarMin,VarMax,nVar,CostFunction)
%% DE Parameters

% nVar= Number of Decision Variables
VarSize=[1 nVar];   % Decision Variables Matrix Size

F=0.3;
CR=0.8;        % Crossover Probability


l=0; 
f=0;
%% Initialization
empty_individual.Position=[];
empty_individual.Cost=[];
BestSol.Cost=inf;
pop=repmat(empty_individual,nPop,1);
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
        S=BestSol.Cost;
        P=BestSol.Position;
    end
    
    f=f+1;
    if f>Max_FEs
        return;
    end
    if mod(f,nPop)==0
        l=l+1;
        Convergence_curve(l)=BestSol.Cost;
    end
end
% BestCost=zeros(MaxIt,1);
%% DE Main Loop
while 1
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
        y=pop(a).Position+F.*(pop(b).Position-pop(c).Position);
        y = max(y, VarMin);
		y = min(y, VarMax);
		
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=CR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        NewSol.Position=z;
        NewSol.Cost=CostFunction(NewSol.Position);
        
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
               S=BestSol.Cost;
               P=BestSol.Position;
            end
        end
        f=f+1;
        if f>Max_FEs
            return;
        end
        if mod(f,nPop)==0
            l=l+1;
            Convergence_curve(l)=BestSol.Cost;
        end
        
    end
    
    % Update Best Cost
    % BestCost(it)=BestSol.Cost;
   
end
end


