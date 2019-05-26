%%%%% Homework 4 Part a (i) %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/21/2019 %%%%%%%

nstepp=14;
nsteps=2^nstepp;
ntjmax=7;
tmax=1000;
variance=tmax/nsteps;
amp=0;
field=0;
g0=1;
K = (6435*pi)/(16384);
Delta_W_n =zeros(nstepp,nsteps); 
realization=randn(1,nsteps); 
x=zeros(ntjmax,nsteps);
Delta_W_n(1,:) = realization(1,:)*sqrt(variance); % for W_n_0

for k=2:7
    Delta_W_n(k,1:2^(nstepp-k+1)) = Delta_W_n(k-1,1:2:length(Delta_W_n(2,1:2^(nstepp-(k-1)+1)))-1)+ Delta_W_n(k-1,2:2:length(Delta_W_n(k-1,1:2^(nstepp-(k-1)+1))));
end

Standard_deviation_vec = zeros(ntjmax,1);
mean_vec = zeros(ntjmax,1);

for ntj = ntjmax:-1:1
        ntfactor=2^(ntj-1);
        nt=nsteps/ntfactor;
        dt(ntj)=tmax/nt;
        for i=1:nt
        x(ntj,i+1)=x(ntj,i)+dt(ntj)*(((cos(x(ntj,i))).^16)/K - 1/(2*pi))+g0*Delta_W_n(ntj,i);
        end
        mean_vec(ntj,1) = mean(x(ntj,:));
        Standard_deviation_vec(ntj,1) = std(x(ntj,:));
        plot(0:dt(ntj):tmax,x(ntj,1:nt+1));
        xlabel('time');
        ylabel('x(t)');
        hold on
end

legend(num2str(dt(ntjmax)),num2str(dt(ntjmax-1)),num2str(dt(ntjmax-2)),num2str(dt(ntjmax-3)),num2str(dt(ntjmax-4)),num2str(dt(ntjmax-5)));
