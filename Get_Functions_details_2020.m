%% CEC 2020
% Dim=10/20

function [lb,ub,dim,fobj] = Get_Functions_details_2020(F,dim)

switch F
    case 'F1'

        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',1); 
        
    case 'F2'

        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',2); 
        
    case 'F3'

        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',3); 
        
        
    case 'F4'
        
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',4); 
        
    case 'F5'
        
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',5); 
        
        
    case 'F6'
        
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',6); 
        
    case 'F7'
        
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',7); 
        
    case 'F8'
        
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',8); 
        
    case 'F9'
        
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',9); 
        
    case 'F10'
        
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec20_func(x',10);     
  
        
end

end
