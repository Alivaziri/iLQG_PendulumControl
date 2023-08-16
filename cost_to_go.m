function cost=cost_to_go(x,u)
cost = 0;

for k = 1 : 400  
    
        if k < 400

        cost = cost + 0.01*cost_calculator1(u(:,k));
        
        else
            
        cost = cost + cost_calculator2(x(:,400));
        
        end 
        
end
