%%%%% Homework 4 Part a(ii) %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/21/2019 %%%%%%%

tic
M_trials = 100; % define the trials numbers
nstepp=14;  % define the factor of time step
nsteps=2^nstepp; % define the total size of the time step
tmax=2000;  % define the max time 
variance=tmax/nsteps; % define the delta t
nsteps_100 = round(nsteps * (100/2000)); % define the T=100
g0=0.2; % sigma =0.2 in this question
realization=randn(M_trials,nsteps); % define the W_0
x=zeros(M_trials,nsteps);  % define the matrix of X(t), which have M trials rows and nsteps time step columns 
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

%%%%%%%%%%%%%%%%% calculate the V(t) %%%%%%%%%%%

V= (-1/(12870*pi))*(5720*sin(2.*x)+2002*sin(4.*x)+728*sin(6.*x)+(455/2)*sin(8.*x) ...
    + 56*sin(10.*x)+10*sin(12.*x)+(8/7)*sin(14.*x)+(1/16)*sin(16.*x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(x(1,:),V(1,:))
ylabel('$V(X(t))$','Interpreter','latex','FontSize',13)
xlabel('$X(t)$','Interpreter','latex','FontSize',13)
title('the motion of the particle in the potential V(x)')

mean_x_T = mean(x(:,end))/tmax;

figure;
h = histogram(V(:,nsteps_100:end));
title('histogram of the potential energy V (X)')
figure;
semilogy(h.BinEdges(1:end-1),h.Values)
title('linlog scale plot of the histogram of the potential energy V (X)')

toc