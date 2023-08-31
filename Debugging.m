   clear all
   clc
S0 = 1.0;   
   
B_up4 = S0 *0.85;
B_down= S0 *0.65;
alpha = 0.0377; %% 3.77% return
N_x_node=1000 ;
N_t_node=1000 ;
x_ini = 0;
x_term = 3;
t_ini = 0;
t_term = 1;

    h_x = (x_term-x_ini)/(N_x_node - 1) ; %% This is the distance between the nodes. We assume equi-distance.
    h_t = (t_term-t_ini)/(N_t_node - 1) ; %% This is the distance between the nodes. We assume equi-distance.

    x = linspace(x_ini + h_x,x_term,N_x_node-1); %% Set the subintervals(partition points) of open ball (0 + h(j),1 - h(j)) according to the number of nodes.
    t = linspace(t_ini, t_term, N_t_node); 
    
    Init_cond_KI = zeros(N_x_node-1,1);
    Init_cond_NOKI= zeros(N_x_node-1,1);
    
    for i = 1:N_x_node-1
        
        if x(i) >= B_down
        
            Init_cond_NOKI(i) = alpha * S0;
               
        elseif x(i)< B_down
            
            Init_cond_NOKI(i) = x(i) - S0 ;
            
        end
        
    end
    
     for i = 1:N_x_node-1
        
        if x(i) >= B_up4
        
            Init_cond_KI(i) = alpha * S0;
               
        elseif x(i)< B_up4
            
            Init_cond_KI(i) = x(i) - S0;
            
        end
        
    end
    plot(x,Init_cond_NOKI,'r-','LineWidth',2)
    hold on
    
    plot(x,Init_cond_KI,'b-','LineWidth',2)
    
    grid on