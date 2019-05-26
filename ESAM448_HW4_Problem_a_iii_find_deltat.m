nstepp=17;
nsteps=2^nstepp;
ntjmax=7;
tmax=2000;
variance=tmax/nsteps;
amp=0;
field=0;
% Parameters
% g=g0+g1*x
g0=0.2;
g1=0;
posit = 1: nsteps; % define the index variables
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
        x(ntj,i+1)=x(ntj,i)+dt(ntj)*F(x(ntj,i),amp,field)+(g0+g1*x(ntj,i))*Delta_W_n(ntj,i);
        end
        mean_vec(ntj,1) = mean(x(ntj,:));
        Standard_deviation_vec(ntj,1) = std(x(ntj,:));
        plot(0:dt(ntj):tmax,x(ntj,1:nt+1));
        xlabel('time');
        ylabel('x(t)');
        hold on
end

legend(num2str(dt(ntjmax-1)),num2str(dt(ntjmax-2)),num2str(dt(ntjmax-3)),num2str(dt(ntjmax-4)),num2str(dt(ntjmax-5)));
