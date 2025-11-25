%% GJO algorithm
function [Male_Jackal_score,Male_Jackal_pos,Convergence_curve]=GJO(SearchAgents_no,Max_FEs,lb,ub,dim,fobj)
    Male_Jackal_pos=zeros(1,dim);
    Male_Jackal_score=inf; 
    Female_Jackal_pos=zeros(1,dim);  
    Female_Jackal_score=inf;

    Positions=initalization(SearchAgents_no,dim,ub,lb);
    
    l=0; 
    f=0;

    while 1

        for i=1:size(Positions,1)

            Flag4ub=Positions(i,:)>ub;
            Flag4lb=Positions(i,:)<lb;
            Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

            fitness=fobj(Positions(i,:));

            if fitness<Male_Jackal_score
                Male_Jackal_score=fitness;
                Male_Jackal_pos=Positions(i,:);
            end

            if fitness>Male_Jackal_score && fitness < Female_Jackal_score
                Female_Jackal_score=fitness;
                Female_Jackal_pos=Positions(i,:);
            end

            f=f+1;
            if f>Max_FEs
                return;
            end
            if mod(f,SearchAgents_no)==0
                l=l+1;
                Convergence_curve(l)=Male_Jackal_score;
            end
        end

        E1=1.5*(1-(f/Max_FEs));
        RL=0.05*levy(SearchAgents_no,dim,1.5); % [-2,2]

        for i=1:size(Positions,1)
            for j=1:size(Positions,2)
                E0=2*rand()-1;
                E=E0*E1;

                if abs(E)<1 
                    D_male_jackal=abs((RL(i,j)*Male_Jackal_pos(j)-Positions(i,j))); 
                    Male_Positions(i,j)=Male_Jackal_pos(j)-E*D_male_jackal;

                    D_female_jackal=abs((RL(i,j)*Female_Jackal_pos(j)-Positions(i,j))); 
                    Female_Positions(i,j)=Female_Jackal_pos(j)-E*D_female_jackal;
                else 
                    D_male_jackal=abs((Male_Jackal_pos(j)-RL(i,j)*Positions(i,j))); 
                    Male_Positions(i,j)=Male_Jackal_pos(j)-E*D_male_jackal;

                    D_female_jackal=abs((Female_Jackal_pos(j)-RL(i,j)*Positions(i,j))); 
                    Female_Positions(i,j)=Female_Jackal_pos(j)-E*D_female_jackal;
                end
                Positions(i,j)=(Male_Positions(i,j)+Female_Positions(i,j))/2; 
                
            end
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