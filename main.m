clear

D=30;
SearchAgents_no=30; 
Max_FEs=10000*D; % Maximum number of function evaluations

Function_name='F1';

disp(['Function:',Function_name]);

[lb,ub,dim,fobj]=Get_Functions_details_2017(Function_name,D);

runn=30;
cost=zeros(runn,1);
run_time=zeros(runn,1);

pos=zeros(runn,dim);

disp(['Max_FEs=',num2str(Max_FEs),',PopulationSize=',num2str(SearchAgents_no),',Runs=',num2str(runn)]);

for i=1:runn
    tic;
    [Male_Jackal_score,Male_Jackal_pos,GJO_cg_curve]=GJO(SearchAgents_no,Max_FEs,lb,ub,dim,fobj);
    run_time(i,:)=toc;
    pos(i,:)=Male_Jackal_pos;
    cost(i,:)=Male_Jackal_score;
    curve(i,:)=GJO_cg_curve;
end
GJO_cg_curve=mean(curve);
disp(['GJO','  ',num2str(min(cost)),'  ',num2str(mean(cost)),'  ',num2str(std(cost)),'  ',num2str(mean(run_time))]);

for i=1:runn
    tic;
    [Male_Jackal_score,Male_Jackal_pos,PSO_cg_curve]=PSO(SearchAgents_no,Max_FEs,lb,ub,dim,fobj);
    run_time(i,:)=toc;
    pos(i,:)=Male_Jackal_pos;
    cost(i,:)=Male_Jackal_score;
    curve(i,:)=PSO_cg_curve;
end
PSO_cg_curve=mean(curve);
disp(['PSO','  ',num2str(min(cost)),'  ',num2str(mean(cost)),'  ',num2str(std(cost)),'  ',num2str(mean(run_time))]);

for i=1:runn
    tic;
    [Male_Jackal_score,Male_Jackal_pos,DE_cg_curve]=DE(SearchAgents_no,Max_FEs,lb,ub,dim,fobj);
    run_time(i,:)=toc;
    pos(i,:)=Male_Jackal_pos;
    cost(i,:)=Male_Jackal_score;
    curve(i,:)=DE_cg_curve;
end
DE_cg_curve=mean(curve);
disp(['DE','  ',num2str(min(cost)),'  ',num2str(mean(cost)),'  ',num2str(std(cost)),'  ',num2str(mean(run_time))]);

for i=1:runn
    tic;
    [Male_Jackal_score,Male_Jackal_pos,QLGJO_cg_curve]=QLGJO(SearchAgents_no,Max_FEs,lb,ub,dim,fobj);
    run_time(i,:)=toc;
    pos(i,:)=Male_Jackal_pos;
    cost(i,:)=Male_Jackal_score;
    curve(i,:)=QLGJO_cg_curve;
end
QLGJO_cg_curve=mean(curve);
disp(['QLGJO','  ',num2str(min(cost)),'  ',num2str(mean(cost)),'  ',num2str(std(cost)),'  ',num2str(mean(run_time))]);

for i=1:runn
    tic;
    [Male_Jackal_score,Male_Jackal_pos,SequoiaOA_cg_curve]=SequoiaOA(SearchAgents_no,Max_FEs,lb,ub,dim,fobj);
    run_time(i,:)=toc;
    pos(i,:)=Male_Jackal_pos;
    cost(i,:)=Male_Jackal_score;
    curve(i,:)=SequoiaOA_cg_curve;
end
SequoiaOA_cg_curve=mean(curve);
disp(['SequoiaOA','  ',num2str(min(cost)),'  ',num2str(mean(cost)),'  ',num2str(std(cost)),'  ',num2str(mean(run_time))]);

for i=1:runn
    tic;
    [Male_Jackal_score,Male_Jackal_pos,H5N1_cg_curve]=H5N1(SearchAgents_no,Max_FEs,lb,ub,dim,fobj);
    run_time(i,:)=toc;
    pos(i,:)=Male_Jackal_pos;
    cost(i,:)=Male_Jackal_score;
    curve(i,:)=H5N1_cg_curve;
end
H5N1_cg_curve=mean(curve);
disp(['H5N1','  ',num2str(min(cost)),'  ',num2str(mean(cost)),'  ',num2str(std(cost)),'  ',num2str(mean(run_time))]);

for i=1:runn
    tic;
    [Male_Jackal_score,Male_Jackal_pos,BKA_cg_curve]=BKA(SearchAgents_no,Max_FEs,lb,ub,dim,fobj);
    run_time(i,:)=toc;
    pos(i,:)=Male_Jackal_pos;
    cost(i,:)=Male_Jackal_score;
    curve(i,:)=BKA_cg_curve;
end
BKA_cg_curve=mean(curve);
disp(['BKA','  ',num2str(min(cost)),'  ',num2str(mean(cost)),'  ',num2str(std(cost)),'  ',num2str(mean(run_time))]);

for i=1:runn
    tic;
    [Male_Jackal_score,Male_Jackal_pos,QTGJO_cg_curve]=QTGJO(SearchAgents_no,Max_FEs,lb,ub,dim,fobj);
    run_time(i,:)=toc;
    pos(i,:)=Male_Jackal_pos;
    cost(i,:)=Male_Jackal_score;
    curve(i,:)=QTGJO_cg_curve;
end
QTGJO_cg_curve=mean(curve);
disp(['QTGJO','  ',num2str(min(cost)),'  ',num2str(mean(cost)),'  ',num2str(std(cost)),'  ',num2str(mean(run_time))]);


%% Figure

figure1=figure('Color',[1 1 1]); 
CNT=20;
Max_iteration=Max_FEs/SearchAgents_no; 
k=round(linspace(1,Max_iteration,CNT)); 
iter=SearchAgents_no:SearchAgents_no:Max_FEs;

hold on
plot(iter,GJO_cg_curve,'Color',"#19161a",'lineWidth',1.5,'Marker','o','MarkerSize',8,'MarkerIndices', k);
hold on
plot(iter,PSO_cg_curve,'Color',"#FF6F71",'lineWidth',1.5,'Marker','square','MarkerSize',8,'MarkerIndices', k);
hold on
plot(iter,DE_cg_curve,'Color',"#731700",'lineWidth',1.5,'Marker','x','MarkerSize',8,'MarkerIndices', k);
hold on
plot(iter,QLGJO_cg_curve,'Color',"#1C7DAA",'lineWidth',1.5,'Marker','*','MarkerSize',8,'MarkerIndices', k);
hold on
plot(iter,SequoiaOA_cg_curve,'Color',"#DE6703",'lineWidth',1.5,'Marker','x','MarkerSize',8,'MarkerIndices', k);
hold on
plot(iter,H5N1_cg_curve,'Color',"#ED2F93",'lineWidth',1.5,'Marker','pentagram','MarkerSize',8,'MarkerIndices', k);
hold on
plot(iter,BKA_cg_curve,'Color',"#33A0CE",'lineWidth',1.5,'Marker','^','MarkerSize',8,'MarkerIndices', k);
hold on
plot(iter,QTGJO_cg_curve,'Color',"#6EAC25",'lineWidth',2,'Marker','hexagram','MarkerSize',10,'MarkerIndices', k);

xlabel('FEs');
ylabel('Fitness');
box on

legend('GJO','PSO','DE','QLGJO','SequoiaOA','H5N1','BKA','QTGJO','Location', 'best')
set(gca,'Yscale','log');

