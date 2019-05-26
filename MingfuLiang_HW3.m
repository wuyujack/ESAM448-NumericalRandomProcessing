%%%%% Homework 3(a) %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/05/2019 %%%%%%%

function [t_exit]=MingfuLiang_HW3(n_trial,L)
tic

%L = 20;
%n_trial = 1000;
%L_2 = L*L;

%%%%%%%%%%%% initialize the parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%

N_trials = n_trial;
r_inject = 0.5;
r_remove = 1.1;
N_S = zeros(1000,400); % define the vector to store the number of susceptible individual for each iteration
N_I = zeros(1000,400); % define the vector to store the number of infected individual for each iteration
N_R = zeros(1000,400); % define the vector to store the number of removed individual for each iteration
t_exit_mat = zeros(1000,1); % define the vector to store the exit time for each iteration
posit = 1: L; % define the index variables 
up_shift = circshift(posit,1); % shift the variables up one unit
down_shift = circshift(posit,-1); % shift the variables down one unit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i =1: N_trials
    
    % Starting with an initial condition in which only a single individual in the center of the
    % system (i = 10; j = 10) is infected (L = 20)
    % here denote 1 is susceptible, 2 is infected, 3 is removed
    
    area_mat = ones(L, L); % define the matrix as the system
    area_mat(10,10) = 2; % define the (10,10) in the system to be susceptiable
    t_counter =1;
    t_exit = 0;
    
    % initialize the NS, NI, NR
    NS = 399;
    NI = 1;
    NR = 0;
    
    %%%%%%%%%%%%% uncomment the code below so that you can watch a moive : )
    
    %figure;
    %h11=imagesc(area_mat);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        while NI~=0
                       % define the infected neighbor matrix
                       Infected_Neigh=area_mat==2;
                       
                       % define the transition rate matrix
                       transit_rates = r_inject*( 0 + Infected_Neigh(up_shift, posit) ...
                                                                          +Infected_Neigh(down_shift, posit) ...
                                                                             +Infected_Neigh(posit, up_shift) ...
                                                                          +Infected_Neigh(posit, down_shift)  ...
                                                            +Infected_Neigh(up_shift, down_shift) ...
                                                         +Infected_Neigh(down_shift, down_shift) ...
                                                                   +Infected_Neigh(up_shift, up_shift)  ...
                                                             +Infected_Neigh(down_shift, up_shift));
                       transit_rates(area_mat==2) =r_remove;
                       transit_rates(area_mat==3) =0;

            % define the sigma c_km, which is the sum of c_km for m =1, m~=j
            sigma_c_km= sum(transit_rates(:));
            
            % define the tau
            tau = -log(rand())/sigma_c_km;
            
            %define the normalized transition rate matrix
            normedrates = transit_rates/sigma_c_km;
            normedrates_vec = normedrates(:);
            
            % find the smallest j satisfies the (sigma ckm from 1 to j) <=
            % U2 * (sigma ckm from m =1 and m~=k)
            i_select=find(cumsum(normedrates_vec)>=rand(), 1);
            
            area_mat_vec = area_mat(:);% stack the matrix to vector
            
            % update the NS, NI, NR
            if area_mat_vec(i_select) ==2
                area_mat_vec(i_select) =3;
                NR = NR +1;
                NI = NI - 1;
            elseif area_mat_vec(i_select) ==1
                area_mat_vec(i_select) = 2;
                NI = NI + 1;
                NS = NS - 1;
            end
            
           area_mat = reshape(area_mat_vec,L,L);% reshape vector into matrix back
        
        % update t_exit and store the stats at this iteration
        t_exit = t_exit  + tau;   
        N_S(i,  t_counter) = NS;
        N_I(i,  t_counter) = NI;
        N_R(i,  t_counter) = NR;
        t_counter = t_counter +1;
        end
        
    t_exit_mat(i,1)=t_exit;
    
    %%%%%%%%%%%%% uncomment the code here so that you can watch a moive : )
    
    %set(h11,'CData',area_mat)
    %pause(0.01);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
figure;
histogram(t_exit_mat,20);
title1 = 'Histogram of mean of t exit';
title(title1);
mean_t_exit = mean(t_exit_mat);
text = ['mean of t_exit is ', num2str(mean_t_exit),' \n'];
fprintf(text);
t_exit = t_exit_mat;
toc
end
