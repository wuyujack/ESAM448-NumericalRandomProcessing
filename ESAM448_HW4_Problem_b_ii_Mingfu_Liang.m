%%%%% Homework 4 Part b (ii) %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/21/2019 %%%%%%%

M_trials =1000; % define the trials numbers
tau = 10;
sigma = 0.5;
nstepp=13; % define the factor of time step
nsteps=2^nstepp; % define the total size of the time step
tmax=1000; % define the max time 
variance=tmax/nsteps; % define the delta t
amp=0;
field=0;
g0=1;
K = (6435*pi)/(16384);

posit = 1: nsteps; % define the index variables
realization=randn(M_trials,nsteps); % define the W_0
x=zeros(M_trials,nsteps); % define the matrix of X(t), which have M trials rows and nsteps time step columns 
Delta_W_n = realization*sqrt(variance); % for W_n_0
G = zeros(M_trials,nsteps); % define the matrix of the G(t), which have the same size of X(t)

ntj = 1;
ntfactor=2^(ntj-1);
nt=nsteps/ntfactor;
dt(ntj)=tmax/nt;

%%%%%%%%%%%%% compute all of the samples simultaneously
%%%%%%%%%%%%% using a single time-stepping loop.

for i=1:nt
        G(:,i+1)= G(:,i)+(-1/tau).*G(:,i).*dt(ntj)+(1/sqrt(tau)).*sigma.*Delta_W_n(:,i);
        x(:,i+1)=x(:,i)+dt(ntj).*(((cos(x(:,i))).^16)/K - 1/(2*pi))+g0.*G(:,i).*dt(ntj);
end

%%%%%%%%%%%%%%%%%%%%%%%%% vectorizing calculate V(X(t)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V= (-1/(12870*pi))*(5720*sin(2.*x)+2002*sin(4.*x)+728*sin(6.*x)+(455/2)*sin(8.*x) ...
    + 56*sin(10.*x)+10*sin(12.*x)+(8/7)*sin(14.*x)+(1/16)*sin(16.*x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(x(1,:),V(1,:))
ylabel('$V(X(t))$','Interpreter','latex','FontSize',13)
xlabel('$X(t)$','Interpreter','latex','FontSize',13)
title('the motion of the particle in the potential V(x)')

%%%%%%% plot the autocorrelation averaging over large trials %%%%%%%%%%%%%%

mean_x = mean(x);
G_cov=zeros(M_trials,2*100+1);
for i =1:M_trials
G_cov(i,:) = xcov(G(i,:),100);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(mean(G_cov))
title(['T = ',num2str(tmax),', M trials =', num2str(M_trials), ' , Tau = ', num2str(tau),' sigma = ',num2str(sigma)])