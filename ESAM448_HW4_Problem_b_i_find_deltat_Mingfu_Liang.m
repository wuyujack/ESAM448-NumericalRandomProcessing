%%%%% Homework 4 Part b(i) %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/21/2019 %%%%%%%

M_trials =10000; % define the trials numbers
tau = 10;
sigma = 0.5; 
nstepp=10; % define the factor of time step
nsteps=2^nstepp; % define the total size of the time step
ntjmax=10; % define how many different delta t I would want to evaluate
tmax=2; % define the max time 
g0=1;
variance=tmax/nsteps; % define the delta t
posit = 1: nsteps; % define the index variables
Delta_W_n =zeros(nstepp,nsteps,M_trials);  
realization=randn(1,nsteps,M_trials); % define the W_0

 % define the 3 dimension matrix of X(t), which have M trials 2 dimension matrix, in each two dimension matrix
 % it has ntjmax row denote different each delta t and nsteps column for each delta t
x=zeros(ntjmax,nsteps,M_trials);
                                                      
G = zeros(ntjmax,nsteps,M_trials); % define the matrix of the G(t), which have the same size of X(t)
Delta_W_n(1,:,:) = realization(1,:,:)*sqrt(variance); % for W_n_0
K = (6435*pi)/(16384);

%%%%%%%%%%%%%%%%%%%% calculate W_n %%%%%%%%%%%%%%%%%%%%%%%%%%

for k=2:ntjmax
    Delta_W_n(k,1:2^(nstepp-k+1),:) = Delta_W_n(k-1,1:2:length(Delta_W_n(2,1:2^(nstepp-(k-1)+1)))-1,:)+ Delta_W_n(k-1,2:2:length(Delta_W_n(k-1,1:2^(nstepp-(k-1)+1))),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_M_trials_X_t = zeros(ntjmax,nsteps+1); % define the mean of X(t) over M trials

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% loop for all different delta t and find out which one would have
%%%%%%% smallest weak error

for ntj = ntjmax:-1:1
        ntfactor=2^(ntj-1);
        nt=nsteps/ntfactor;
        dt(ntj)=tmax/nt;
        for i=1:nt
        G(ntj,i+1,:)= G(ntj,i,:) + (-1/tau).*G(ntj,i,:).*dt(ntj)+(1/sqrt(tau)).*sigma.*Delta_W_n(ntj,i,:);
        x(ntj,i+1,:)=x(ntj,i,:)+dt(ntj).*(((cos(x(ntj,i,:))).^16)/K - 1/(2*pi))+g0.*G(ntj,i,:).*dt(ntj);
        end
        M_trials_X_t = reshape(x(ntj,1:nt+1,:),[nt+1,M_trials]).';
        mean_M_trials_X_t(ntj,1:nt+1) = mean(M_trials_X_t(:,1:nt+1));
        mean_M_trials_X_T(ntj,1) = mean_M_trials_X_t(ntj,nt+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weak_error_mean_shift =abs(mean_M_trials_X_T(2:end,1) -mean_M_trials_X_T(1:end-1,1));

figure();
loglog(dt,dt);
hold on
loglog(dt(2:nstepp),weak_error_mean_shift)
ylabel('$ weak \quad diff. \quad between \quad succ. \quad approx.$','Interpreter','latex','FontSize',13)
xlabel('$\Delta {t}$','Interpreter','latex','FontSize',13)
legend('delta t verse delta t','weak error verse delta t')
title(['T = ', num2str(tmax), ', tau =', num2str(tau), ', sigma = ', num2str(sigma), ', M trials = ', num2str(M_trials)])
hold off

