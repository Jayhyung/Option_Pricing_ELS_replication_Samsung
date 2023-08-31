clear all
clc


%%%%%%%%%%%%%%%%% Calculating ELS using FDM Method %%%%%%%%%%%%%%%%%%%%%%%%
% Before doing iterations, we set related parameters, exercise price(E),
% risk-free rate(r),volatility(sigma)

r = 0.03;
sigma=0.25;
S0 = 0.9;
B_up4 = S0 *0.85;
B_down= S0 *0.65;
B_up = [1 1 1 0.95 0.9].*S0 ;
alpha = 0.0377*2; %% 3.77% return
N_x_node=751;
N_t_node= 751 ;
x_ini = 0;
x_term = 3;
t_ini = 0;
t_term = 3;
x = linspace(x_ini,x_term,N_x_node);
t = linspace(t_ini,t_term,N_t_node);
Monitoring_Period = [0.5, 1.0, 1.5, 2.0, 2.5]; %% We do early redemption check 5 times.
                                           % In fact, final check is
                                           % already described in the
                                           % initial condition when solving
                                           % pde.
  

[h_x,h_t,true_x,true_t, u_KI, u_NOKI,Indicator_KI,Indicator_NOKI,Init_cond_KI,Init_cond_NOKI,u_aux] = Black_Scholes_Implicit_B(r,sigma,N_x_node,N_t_node, Monitoring_Period, x_ini,x_term,t_ini,t_term,S0,B_down,B_up,B_up4,alpha);

delta_ELS = zeros(N_x_node-2, N_t_node);

for i = 1:N_x_node-2
   
    for j = 1:N_t_node
        if i ==1 || N_x_node-2
        delta_ELS(i,j)=(u_NOKI(i+1,j)-u_NOKI(i,j))/h_x;
        else
        delta_ELS(i,j)=(u_NOKI(i+1,j)-u_NOKI(i-1,j))/(2*h_x);    
        end
    end
    
end

      figure(1)
        [t,x]=meshgrid(true_t,true_x); %% Just save t = [0.9,1.0]
   
        mesh(x,t,u_NOKI)
        hold on
        
        xlabel('x')
        ylabel('t')
        zlabel('u(x,t)')
        
        title('Numerical Solution of ELS')
        grid on
     
     figure(2)
     [t,x]=meshgrid(true_t,true_x(1:end-1));
        mesh(x,t,delta_ELS)
        hold on
        
        xlabel('x')
        ylabel('t')
        zlabel('Delta')
        title('Numerical Solution of ELS Delta')
        grid on
        zlim([-1 5])
        
        u_NOKI(find(true_x==1),end)