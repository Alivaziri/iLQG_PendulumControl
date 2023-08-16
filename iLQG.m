clc;close all;clear all;
%choose initial conditions (u,x0)
u0=1;
x0=[0;0];
x1=zeros(2,400);
u=ones(1,399);   
A=zeros(2,2,399);
Max_Iteration=100;

% First iteration, forward simulation using dynamic equation 1
for k = 1 : 399      
        x1(:,k+1) =x1(:,k)+0.01* frw(x1(:,k),u(:,k));
        x_new(:,k+1)=x1(:,k+1);
        cost = cost_to_go(x1 , u); %%% cost-to-go function is used for calculating the cost using position matrix and control input vector at the first iteration 
end

for i = 1 : Max_Iteration
    %%% calculating s0,s,S in the final time step in each iteration to go
    %%% backward for calculating g, G, H and other s0, s, S for the whole
    %%% time. There can be find a final cost using final position with which s0, s, S are calculated
    
    ii=i; %%counter of iterations
    s0(1,400) = 100*(1+cos(x_new(1,400)))^2 + 10*(x_new(2,400))^2;
    s(: ,:, 400) = [-200*sin(x_new(1,400))*(1 + cos(x_new(1,400))) ; 20*x_new(2,400)];
    S(: ,:, 400) = [-200*(cos(x_new(1,400)) + cos(2*x_new(1,400))),0 ;0 ,20];
    
    for k = 399 : -1 : 1  
        %%% Calculation of matrics for the nonlinear discritized dynamics
        B = [0;0.01];
        R=0.02; 
        q0(k) = 0.01*sum(u(k).^2);
        r(:,k) = 0.02*u(k);
        g = r(:,k) + B'*s(:,k+1);
        G =  B'*S(:,:,k+1)*A(:,:,k);
        H = R + B'*S(:,:,k+1)*B;
 
        l(:,k) = -inv(H)*g;
        L(:,:,k) = -inv(H)*G;
        
        if abs(x0(1))> 0
            l(:,k) = min(max(l(:,k)+u(:,k),-2),2) - u(:,k);
            L((l(:,k)+u(:,k)==-2)|(l(:,k)+u(:,k)==2),:,k) = 0;
        end

        A(:,:,k) = [ 1 , 0.01 ;-0.04 * cos(x_new(1 , k)) , 1];
        S(:,:,k) =   A(:,:,k)'*S(:,:,k+1)*A(:,:,k) + L(:,:,k)'*H*L(:,:,k) + L(:,:,k)'*G + G'*L(:,:,k);
        s(:,k) =  A(:,:,k)'*s(:,k+1) + L(:,:,k)'*H*l(:,k) + L(:,:,k)'*g + G'*l(:,k);
        s0(k) = q0(k) + s0(k+1) + l(:,k)'*H*l(:,k)/2 + l(:,k)'*g;
    end
    
    x1(:,1)=x0;
    delta_x = zeros(2,1);
    
    for k = 1:399
        %%% affine control law and its cost using the linearized dynamics,
        %%% I checked this step to calculate cost and position matrix (x)
        %%% using nonlinear dynamics, and the method didn't work! we just
        %%% can use this step along with linearized dynamic system in order
        %%% to get the same results as that of the paper
        delta_u = l(:,k) + L(:,:,k)*delta_x;   
        delta_u = min(max(delta_u+u(:,k),-2),+2) - u(:,k);
        delta_x = A(:,:,k)*delta_x + B*delta_u;  
        u_tilt(:,k) = u(:,k) + delta_u;
        x1(:,k+1) =x1(:,k)+0.01* frw(x1(:,k),u_tilt(:,k));
        cost_new = cost_to_go(x1 , u);
    end

    if cost_new < cost
        u = u_tilt;
        x_new = x1;
    end
    
    if  abs(cost_new - cost)< 0.00001 * cost
        cost = cost_new;
        break;
    end
    cost = cost_new;
end

fprintf('Iteration = %d; \n', ii);
fprintf('OptimalCost =%f \n', cost_new);
figure
plot(x_new(1,:),x_new(2,:));
xlabel('angular position'); ylabel('angular velocity'); grid on;
figure
plot(1:1:399,u);
xlabel('time'); ylabel('control input'); grid on;
