function [h_x,h_t,x,t, u_KI, u_NOKI,Indicator_KI,Indicator_NOKI,Init_cond_KI,Init_cond_NOKI,u_aux] = Black_Scholes_Implicit_B(r,sigma,N_x_node,N_t_node, Monitoring_Period, x_ini,x_term,t_ini,t_term,S0,B_down,B_up,B_up4,alpha)
%% Let's solve the Black_Scholes equation numerically using Implicit scheme.


%% Explanation of Input variables
% "r" is risk-free rate
% "sigma" is volatility
% "E" is exercise price of financial option.
% "Call_Put_Type" is option type indicator.
% "N_x_node" is the # of x nodes.
% "N_t_node" is the # of t nodes.
% "B_up" is the boundary value of barrier option
% "Monitoring_Period" is the monitoring period of barrier option
% "x_ini" is the initial value of x space.
% "x_term" is the terminal value of x space.
% "t_ini" is the initial value of t space.
% "t_term" is the terminal value of t space.


%% Explanation of Output variables
% "h_x" is the distance value of x nodes.
% "h_t" is the distance value of t nodes.
% "x" is the x interval.
% "t" is the t interval



%% Explanation of Local Variables 
% "u_NOKI" is the solution of "no knocked-in" case. 
% "u_KI"   is the solution of "knocked-in" case.
% "u_aux"  is the auxilary value for calculating solution.
% "Indicator_NOKI" is the indicator function of no knocked-in case(i.e. when x in "NOKI" area, this value is 1, otherwise 0)
% "Indicator_KI" is the indicator function of knocked-in case(i.e. when x in "KI" area, this value is 1, otherwise 0)
% "Init_cond_KI" is the initial condition for "Knocked-in" case.
% "Init_cond_NOKI" is the initial condition for "no Knocked-in case"

           
   %% Define simulated solution vectors,Indicator functions, nod_distance(x and t), and discretize x & t.
    u_NOKI = zeros(N_x_node-1,N_t_node); %% Set the vector for simulated solution attained by FDM.
    u_KI = zeros(N_x_node-1,N_t_node);
    u_aux = zeros(N_x_node-1,N_t_node);
    
    u_n_NOKI = zeros(N_x_node-1,1); %% We use these vectors for saving solution temporarilly.
    u_n_KI   = zeros(N_x_node-1,1);
    u_n_aux = zeros(N_x_node-1,1);
    
    
    Indicator_KI = zeros(N_x_node-1,1);
    Indicator_NOKI = zeros(N_x_node-1,1);
    h_x = (x_term-x_ini)/(N_x_node - 1) ; %% This is the distance between the nodes. We assume equi-distance.
    h_t = (t_term-t_ini)/(N_t_node - 1) ; %% This is the distance between the nodes. We assume equi-distance.
    x = linspace(x_ini + h_x,x_term,N_x_node-1); %% Set the subintervals(partition points) of open ball (0 + h(j),1 - h(j)) according to the number of nodes.
    t = linspace(t_ini, t_term, N_t_node); 
    

    %%%%%%%% Array Indicator functions%%%%%%%%%%%%%%%

    for i = 1:N_x_node-1
   
        if x(i) >=B_down
            Indicator_NOKI(i) = 1; 
        
        elseif x(i) <B_down
            Indicator_KI(i) = 1;
        end
    
    end
    
    

  
    epsil = 10^(-8);
    
    
    %% Define the initial conditions of the targetted ELS(Line 54 ~ 90).
    
    Init_cond_KI = zeros(N_x_node-1,1);
    Init_cond_NOKI= zeros(N_x_node-1,1);
    
    for i = 1:N_x_node-1
        
        if x(i) >= B_down
        
            Init_cond_NOKI(i) = (1+3*alpha) * S0;
               
        elseif x(i)< B_down
            
            Init_cond_NOKI(i) = x(i) ;
            
        end
        
    end
    
     for i = 1:N_x_node-1
        
        if x(i) >= B_up4
        
            Init_cond_KI(i) = (1+3* alpha)* S0;
               
        elseif x(i)< B_up4
            
            Init_cond_KI(i) = x(i);
            
        end
        
    end
    


    %% Initiate iteration for solving Black-Scholes equation.
    
    %% Now, to calculate the solution of each node by using Thomas Algorithm, we set the related elements.
     % Note that we constuct different setting order(different from heat equation setting) of the related elements diag, lower, upper, and force 
     % because of the zero-gamma condition of Black-Scholes equation. 
    
      
            
    % Construct Upper and Lower Diagonal Element.    
    
     upper_diag = zeros(N_x_node-2,1);
     lower_diag = zeros(N_x_node-2,1);
   
        for i= 1:N_x_node-2         
            upper_diag(i) = -0.5*((sigma * x(i)/h_x)^2)*h_t - (0.5*r*h_t*x(i)/(h_x)); %% We denote this term as gamma(i)
        end 
        
        for i= 1:N_x_node-3     %% We will set the value lower_diag(N_x_node-2) later, because of zero-gamma condition.
      
            lower_diag(i) = -0.5*((sigma * x(i+1)/h_x)^2)*h_t + (0.5*r*h_t*x(i+1)/(h_x)); %% We denote this term as alpha(i)
            
        end
      
     
    % Construct Diagonal Elements.
    
    diag = zeros(N_x_node-1,1);
        for i = 1:N_x_node-1

            if i == N_x_node-1
                diag(i)= 1-upper_diag(i-1)/lower_diag(length(lower_diag)-1); % zero-gamma condition.
            else
                diag(i) = 1 +h_t*r + ((sigma*x(i)/h_x)^2)*h_t;  % we denote this as beta(i) 
            end
        
        end    

      lower_diag(N_x_node-2) = -2-(diag(N_x_node-2)/lower_diag(length(lower_diag)-1)); %% zero-gamma condition

      
%% Now, we initiate time iteration. This iteration proceeds as the following; 
% First, we know the initial values of "NOKI" & "KI" solutions because we
% are given with initial conditions(Init_cond_NOKI & Init_cond_KI),
% therefore, no need to calculate "u_n_NOKI" & "u_n_KI" at first. 
% 1. We first calculate the solution of "u_n+1_KI" using previous solution "u_n_KI" by Thomas Algorithm.
% 2. Calculate "u_aux" from the "u_n_NOKI" by Thomas Algorithm. 
% 3. Mix the "u_n+1_KI" & "u_aux", implying that we do the following calculation;
% u_n+1_NOKI(xi) = u_n+1_KI(xi)*Indicator_KI + u_aux(xi)*Indicator_NOKI
% 
% 4. Do the early redemption check both on u_n+1_KI & u_n+1_NOKI. 
% 5. Save the values "u_n+1_KI" & "u_n+1_NOKI" as "u_KI(:,n+1)" & "u_NOKI(:,n+1)" respectively
% Now, Use the "u_n+1_NOKI" & "u_n+1_KI" as "u_n_NOKI" & "u_n_KI" and repeat the algorithm again.
    %% Algorithm process 1

             for i = 1 : N_x_node-1
             u_n_KI(i) = Init_cond_KI(i);
             end
             u_KI(:,1) =u_n_KI; 
             
             for i = 1 : N_x_node-1
             u_n_NOKI(i) = Init_cond_NOKI(i);
             end
             u_NOKI(:,1) =u_n_NOKI; 
             
    %% Algorithm process 2         
for n = 1 : N_t_node-1   

        
     % Construct Force
     
     force_NOKI = zeros(N_x_node-1,1);
     force_KI   = zeros(N_x_node-1,1);
     
     for i = 1:N_x_node-1
         if     i ==1
                force_NOKI(i) = u_n_NOKI(i) ;
                force_KI(i) = u_n_KI(i);
         elseif i == N_x_node-1 
                force_NOKI(i) =  -u_n_NOKI(N_x_node-2)/lower_diag(length(lower_diag)-1) ;
                force_KI(i) =  -u_n_KI(N_x_node-2)/lower_diag(length(lower_diag)-1) ;
                
         else
                force_NOKI(i) = u_n_NOKI(i) ;
                force_KI(i) = u_n_KI(i);
         end
     end
     
     %nsize
     nsize = N_x_node-1;
     
     % Use the Thomas Algorithm to get the solution of the targetted PDE.
     
     u_n_aux=TriD_Sol(diag, upper_diag , lower_diag, force_NOKI, nsize);
     u_n_KI=TriD_Sol(diag, upper_diag , lower_diag, force_KI, nsize); % This is "u_n+1_KI"

  %%  Algorithm process 3  
  
  for i = 1 : N_x_node-1
     u_n_NOKI(i) = u_n_KI(i)*Indicator_KI(i) + u_n_aux(i)*Indicator_NOKI(i); %% This is "u_n+1_NOKI"     
  end
  
  %% Algorithm process 4
  scaled_time_point = t(n);

  for j = 5:-1:1

     if abs(scaled_time_point - Monitoring_Period(6-j))<= epsil;
     
            for i = 1:N_x_node-1
            
                if x(i) >= B_up(j)
                     u_n_NOKI(i) =S0*(1+j*0.5*alpha);
                     u_n_KI(i) =S0*(1+j*0.5*alpha);
                end
            end 
    end
  end
            
  %% Algorithm process 5
     
    u_KI(:,n+1) =u_n_KI;
    u_NOKI(:,n+1)=u_n_NOKI;
    u_aux(:,n+1)=u_n_aux;

end
end

