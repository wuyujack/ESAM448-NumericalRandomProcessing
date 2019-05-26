%%%%% Homework 4 Part a (i) and (ii) %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/21/2019 %%%%%%%

tic
M_trials = 8; % define the trials numbers
nstepp=12; % define the factor of time step
nsteps=2^nstepp; % define the total size of the time step
tmax=64; % define the max time 
variance=tmax/nsteps;  % define the delta t
g0=1;
realization=randn(M_trials,nsteps); % define the W_0
x=zeros(M_trials,nsteps); % define the matrix of X(t), which have M trials rows and nsteps time step columns 
Delta_W_n = realization*sqrt(variance); % for W_n_0
K = (6435*pi)/(16384); 
ntj = 1; 
ntfactor=2^(ntj-1);
nt=nsteps/ntfactor;
dt(ntj)=tmax/nt;

%%%%%%%%%%%%%%%%% calculate the X(t) %%%%%%%%%%%

for i=1:nt
        x(:,i+1)=x(:,i)+dt(ntj).*(((cos(x(:,i))).^16)/K - 1/(2*pi))+g0.*Delta_W_n(:,i);
end

%%%%%%%%%%%% several plot of X(t) %%%%%%%%%%%%

figure;
plot(0:dt(ntj):tmax,x(1,:));
hold on
plot(0:dt(ntj):tmax,x(2,:));
hold on
plot(0:dt(ntj):tmax,x(3,:));
hold on
plot(0:dt(ntj):tmax,x(4,:));
hold on
plot(0:dt(ntj):tmax,x(5,:));
hold on
plot(0:dt(ntj):tmax,x(6,:));
hold on
plot(0:dt(ntj):tmax,x(7,:));
hold on
plot(0:dt(ntj):tmax,x(8,:));
hold off
ylabel('$x(t)$','Interpreter','latex','FontSize',13)
xlabel('$t$','Interpreter','latex','FontSize',13)
title(['several trajectories of X(t) for sigma =1 and M trials = ',num2str(M_trials)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

average_rate = mean(x(:,nsteps+1))/1000;
Standard_deviation = std(x);
Standard_deviation_T = std(x(:,nsteps+1)/1000);
Standard_error_mean = Standard_deviation_T/sqrt(M_trials);
mean_x = mean(x);
figure;
plot(0:dt(ntj):tmax,mean_x);
title(['mean of X(t) verse time for M trials = ', num2str(M_trials)])
figure;
plot(0:dt(ntj):tmax,Standard_deviation)
title(['standard derivation of X(t) verse time for M trials = ', num2str(M_trials)])
toc

fprintf(['Need ', num2str(M_trials), ' to achieve standard error of mean ', num2str(Standard_error_mean)],'\n')
