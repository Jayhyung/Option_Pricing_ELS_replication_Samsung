clear all
clc

randn('state',100)
%%%%%%%%%%%%%%%%% Calculating ELS using Monte-Carlo Method %%%%%%%%%%%%%%%%%%%%%%%%
% Before doing iterations, we set related parameters, exercise price(E),
% risk-free rate(r),volatility(sigma)

r = 0.03;
sigma=0.25;
S0 = 1.1;
B_up1 = S0 * 1.0;
B_up2 = S0 * 0.95;
B_up3 = S0 * 0.90;
B_up4 = S0 *0.85;
B_down= S0 *0.65;
alpha = 0.0377*2; %% 3.77% return
Dt = 1/250;
check_knock_in = 0;
N_path = 100000;
Monitoring_Period = [125 250 375 500 625 750];
pay_off = 0;
T1 = 1/2; T2 = 1; T3 = 1.5; T4 = 2.0; T5 = 2.5; T6 = 3.0;

T = 3;
N = 3/Dt;

%% Make Scenarios

V = zeros(N_path,1); %% Value of ELS

for i = 1:N_path


S = S0*cumprod(exp((r-0.5*sigma^2)*Dt + sigma*sqrt(Dt)*randn(N,1))); %% Make a scenario for underlying asset dynamics. 
Smin = min(S);
    if Smin <= B_down  %% Check whether the underlying asset price touches the bottom level during investment period, B_down... 
                       %  If it touches the bottom during investment
                       %  period and no early redemption, principal loss could occur if the
                       %  terminal asset price(S(T)) is below B_up4.
                       
        check_knock_in = 1;
    end
        for n = 1 : length(S)
            
           
            
          %%%%%%%%%%%%%%%%%  Check wether early redemption occurs or not(Early Redemption Check is from line 46 to 87). %%%%%%%%%%%%%%%%%
                if n == Monitoring_Period(1)
                    
                    if S(n) >= B_up1
                    pay_off = S0*(1+0.5*alpha)*exp(-r*T1);
                    
                    break
                    end
               
                    
                elseif n == Monitoring_Period(2)
                    
                    if S(n) >= B_up1
                    pay_off = S0*(1+alpha)*exp(-r*T2);
                    
                    break
                    end
                    
                elseif n == Monitoring_Period(3)
                    
                    if S(n) >= B_up1
                    pay_off = S0*(1+1.5*alpha)*exp(-r*T3);
                    
                    break
                    end
                    
                elseif n == Monitoring_Period(4)
                   
                    if S(n) >= B_up2
                    pay_off = S0*(1+2*alpha)*exp(-r*T4);
                    
                    break
                    end
                    
                elseif n == Monitoring_Period(5)
                    
                    if S(n) >= B_up3
                    pay_off = S0*(1+2.5*alpha)*exp(-r*T5);
                    
                    break
                    end
            %%%%%%%%%%%%%%%%%  We completed the early redemption check %%%%%%%%%%%%%%%%%       
   
                elseif n == Monitoring_Period(6)  %% No early Redemption occur. Check wether the knock_in event occurs, 
                                                  %  and accorind to that,
                                                  %  calculate the value of ELS.
                                                  
                    
                   if check_knock_in ==1 % When Knock_in event occurs during investment period,
                      
                     if S(n) >= B_up4
                     pay_off = S0*(1+3.0*alpha)*exp(-r*T6);
                    
                     elseif S(n) < B_up4
                        
                      pay_off = S(n)*exp(-r*T6);
                        
                     end
                    
                  
                   elseif check_knock_in ==0 % When Knock_in event does not occur during investment period,
                      
                  
                       pay_off = S0*(1+3.0*alpha)*exp(-r*T6);
                       
                   end
                
                end
            
            
       
        end
        
        
        V(i) = pay_off;
end
   
B_numeric = mean(V)
B_numeric_std = std(V);
conf = [B_numeric - 1.96*B_numeric_std/sqrt(N_path), B_numeric + 1.96*B_numeric_std/sqrt(N_path)]
