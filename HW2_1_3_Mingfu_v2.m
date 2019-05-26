%%%%%%% author: Mingfu Liang
%%%%%%% date: 02/27/2019
%%%%%%% Problem 1(c) 

%%%%%%% initialize the parameters and the matrix used for storing the
%%%%%%% statistics that we need to measure

N_trials=500;
p=9;
alpha = 3/4;
x_max = 2;
I_exact = 3.562937573;

I_hat_matrix = zeros(p,N_trials);
Sigma_hat_I = zeros(p,N_trials);
Sigma_absoluate_error = zeros(p,N_trials);
Sigma_std_error =  zeros(p,N_trials);
Std_var_I_hat = zeros(p,N_trials);

epsilon_I_hat_matrix = zeros(p,1);
epsilon_m_matrix = zeros(p,1);
average_I_hat=zeros(p,1);
Sigma_hat_I_mean = zeros(p,1);
Std_I_hat=zeros(p,1);
M_range=zeros(p,1);

%%% To see how three quantities average of I_hat, average of Delta_hat and
%%% average of std(I_hat) scale with M, it is equal to see how they sacle
%%% with the p since M = 4^p with p up to 9

for i =1:p
    M = 4^i;
    M_range(i,1) = M;
    for trials =1:N_trials
     %%% use M=4^i sample points X_j for the Monte-Carlo integration of each of the integral
     %%% to estimate I_hat
     
     X_j = rand(1,M);
     Integral_alpha = 4*(x_max^(1-alpha));
     X_hat = (Integral_alpha*(1-alpha)*X_j).^(4);
     g_X_j= Integral_alpha*exp(-X_hat);
     
     I_hat = sum(g_X_j)/M;
     Abs_val_error = abs(I_hat-I_exact);
     abs_error=(g_X_j-I_hat).^2;
     std_error_mean =sqrt(sum(abs_error,2)/(M*(M-1))); % estimate the standard error of the mean
     
     I_hat_matrix(i,trials)=I_hat;
     Sigma_absoluate_error(i,trials)=Abs_val_error;
     Sigma_hat_I(i,trials)=std_error_mean;
    end
    
    %%% calculate average_{\hat{I}}, <\hat{\sigma{I}}>, std(\hat{I}),
    %%% \epsilon_{m} and \epsilon_{\hat{I}}
    
    average_I_hat(i,1) = mean(I_hat_matrix(i,:),2);
    Sigma_hat_I_mean(i,1)=mean(Sigma_hat_I(i,:),2);
    Std_I_hat(i,1)=std(I_hat_matrix(i,:));
    epsilon_I_hat_matrix(i,1)= abs(mean(I_hat_matrix(i,:),2)-I_exact);
    epsilon_m_matrix(i,1) = mean(Sigma_absoluate_error(i,:),2);
   
end

%%%%%%% Evaluate the relationship between M and average_{\hat{I}}, <\hat{\sigma{I}}>, std(\hat{I})
p_I=polyfit(log10(M_range),log10(average_I_hat),1);  
p_Std_I_hat=polyfit(log10(M_range),log10(Std_I_hat),1);  
p_Sigma_hat_I_mean=polyfit(log10(M_range),log10(Sigma_hat_I_mean),1);  

%%%%%%% Visualize the result
figure;
loglog(M_range,average_I_hat,'-s');
hold on
loglog(M_range,epsilon_m_matrix,'-s');
hold on
loglog(M_range,Std_I_hat,'-s');
hold on
loglog(M_range,epsilon_I_hat_matrix,'-s');
hold on
loglog(M_range,Sigma_hat_I_mean,'-s');
legend('Average Ihat','Epsilon m','Std(Ihat)','Epsilon Ihat','Average sigmahat(I)')
xlabel('$M$','Interpreter','latex','FontSize',13)
ylabel('$Statistics$','Interpreter','latex','FontSize',13)
mytitle1 = ['log-log plot of all the statistics of Problem 1(c) when alpha =',num2str(alpha)];
title(mytitle1);

figure;
loglog(M_range,average_I_hat,'-s');
mytitle4 = 'log-log plot of average I hat';
title(mytitle4);
xlabel('$M$','Interpreter','latex','FontSize',13)
ylabel('$<\hat{I}>$','Interpreter','latex','FontSize',13)

figure;
loglog(M_range,Sigma_hat_I_mean,'-s');
mytitle3 = 'log-log plot of the average of the standard error of the mean';
title(mytitle3);
xlabel('$M$','Interpreter','latex','FontSize',13)
ylabel('$<\hat{\sigma}_{I}>$','Interpreter','latex','FontSize',13)

figure;
loglog(M_range,Std_I_hat,'-s');
mytitle2 = 'log-log plot of standard deviation of I hat';
title(mytitle2);
xlabel('$M$','Interpreter','latex','FontSize',13)
ylabel('$std(\hat{I})$','Interpreter','latex','FontSize',13)

