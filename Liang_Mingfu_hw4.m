function X = Liang_Mingfu_hw4(M,T,sigma,tau)

%%%%% Homework 4 Part c %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/21/2019 %%%%%%%

%%%%%%%%%%%%%%%% initialize the parameter %%%%%%%%%%%%%
    M_trials = M; % define the trials numbers
    nstepp=13; % define the factor of time step
    nsteps=2^nstepp; % define the total size of the time step
    tmax=T; % define the max time 
    variance=tmax/nsteps; % define the delta t
    K = (6435*pi)/(16384); % define the parameter K
    g0=1; 

    realization=randn(M_trials,nsteps); % define the W_0
    x=zeros(M_trials,nsteps); % define the matrix of X(t), which have M trials rows and nsteps time step columns 
    Delta_W_n = realization*sqrt(variance); % for W_n_0
    G = zeros(M_trials,nsteps); % define the matrix of the G(t), which have the same size of X(t)
    
    ntj =1;
    dt(ntj)=tmax/nsteps;
  
%%%%%%%%%%%%% compute all of the samples simultaneously
%%%%%%%%%%%%% using a single time-stepping loop.
    for i=1:nsteps
        G(:,i+1)= G(:,i)+(-1/tau).*G(:,i).*dt(ntj)+(1/sqrt(tau)).*sigma.*Delta_W_n(:,i);
        x(:,i+1)=x(:,i)+dt(ntj).*(((cos(x(:,i))).^16)/K - 1/(2*pi))+g0.*G(:,i).*dt(ntj);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%% Extract the X(T) %%%%%%%%%%%%%
    X = x(:,end);

 %%%%%%%%%%%%%%%%%% plot X(t) for one realization %%%%%%%%%%%%%%%%%%
 
%  figure;
%  plot(0:dt(ntj):tmax, x(1,:))
%  ylabel('$X(t)$','Interpreter','latex','FontSize',13)
%  xlabel('$t$','Interpreter','latex','FontSize',13)
 
 %%%%%%%%%%%%%%%%%% plot <X(t)> for all realization %%%%%%%%%%%%%%%%%%
 
%  figure;
%  plot(0:dt(ntj):tmax, mean(x))
%  ylabel('$<X(t)>$','Interpreter','latex','FontSize',13)
%  xlabel('$t$','Interpreter','latex','FontSize',13)

%%%%%%%%%%%%%%%%%% autocorrelation of G(t) %%%%%%%%%%

%     G_cov=zeros(M_trials,2*100+1);
%     for i =1:M_trials
%         G_cov(i,:) = xcov(G(i,:),100);
%     end
% 
%     figure;
%     plot(mean(G_cov))
%     title(['xcov of G whenT = ',num2str(tmax),', M trials =', num2str(M_trials), ' , Tau = ', num2str(tau),' sigma = ',num2str(sigma)])

end