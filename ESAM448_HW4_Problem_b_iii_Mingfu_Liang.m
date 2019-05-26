%%%%% Homework 4 Part b (iii) %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/21/2019 %%%%%%%

tic
M_trials =100; % define the trials numbers
tau_head = 100; % define the max of tau
tau_bottom = 0.1; % define the min of tau
delta_tau = 0.1; % define the change of tau
tau_range = round((tau_head-tau_bottom)/delta_tau); % define the range of tau
sigma = 0.5; % define sigma
nstepp=15; % define the factor of time step
nsteps=2^nstepp; % define the total size of the time step
tmax=5000; % define the max time 
variance=tmax/nsteps; % define the delta t
g0=1; 
K = (6435*pi)/(16384);

realization=randn(M_trials,nsteps); % define the W_0
x=zeros(M_trials,nsteps);  % define the matrix of X(t), which have M trials rows and nsteps time step columns 
Delta_W_n = realization*sqrt(variance); % for W_n_0
G = zeros(M_trials,nsteps); % define the matrix of the G(t), which have the same size of X(t)

ntj = 1; 
ntfactor=2^(ntj-1);
nt=nsteps/ntfactor;
dt(ntj)=tmax/nt;
mean_X_T_vec = zeros(tau_range,1);
k =1;

%%%%% evaluate the mean drift shif over different tau %%%%%%%%%%%%%%%%%%%

for tau = tau_bottom:delta_tau:tau_head
    
    for i=1:nt
        G(:,i+1)= G(:,i)+(-1/tau).*G(:,i).*dt(ntj)+(1/sqrt(tau)).*sigma.*Delta_W_n(:,i);
        x(:,i+1)=x(:,i)+dt(ntj).*(((cos(x(:,i))).^16)/K - 1/(2*pi))+g0.*G(:,i).*dt(ntj);
    end
    
    mean_X_T = mean(x(:,end))/tmax;
    mean_X_T_vec(k,1)=mean_X_T;
    k = k+1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
loglog(tau_bottom:delta_tau:tau_head,mean_X_T_vec);
xlabel('$\tau$','Interpreter','latex','FontSize',14)
ylabel('$<X(T)>/T$','Interpreter','latex','FontSize',14)
toc