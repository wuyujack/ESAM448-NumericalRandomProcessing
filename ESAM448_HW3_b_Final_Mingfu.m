%%%%% Homework 3(b) %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/05/2019 %%%%%%%

%%%%%%%%%%%% initialize the parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
L =100;
L_2 = L*L;
N_trials = 1;
r_inject = 0.5;
N_S = zeros(6,1); % define the vector to store the number of susceptible individual for each iteration
N_I = zeros(6,1); % define the vector to store the number of infected individual for each iteration
N_R = zeros(6,1); % define the vector to store the number of removed individual for each iteration
t_exit_mat = zeros(6,1); % define the vector to store the exit time for each iteration
k=1;
posit = 1: L; % define the index variables 
up_shift = circshift(posit,1); % shift the variables up one unit
down_shift = circshift(posit,-1); % shift the variables down one unit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r_remove = 0:0.3:1.5
    
    % Starting with an initial condition in which only a single individual in the center of the
    % system (i = 10; j = 10) is infected (L = 20)
    % here denote 1 is susceptible, 2 is infected, 3 is removed
    
    area_mat = ones(L, L); % define the matrix as the system
    area_mat(50,47:53) = 2; % define the (50,47:53) in the system to be susceptiable
    t_counter =1;
    t_exit = 0;
    
    %%%%%%%%%%%%% uncomment the code below so that you can watch a moive : )
    
    %figure;
    %h11=imagesc(area_mat);
    %title1 = ['snapshot of r remove = ',num2str(r_remove)];
    %title(title1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initialize the NS, NI, NR
    NS = L_2-7;
    NI = 7;
    NR = 0;
    
        while NI~=0
            if NI == L_2
                break
            end
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
            
            % define the random number and tau 
            U1 = rand();
            tau = -log(U1)/sigma_c_km;
            U2 = rand();
            sigma_ckm = 0;
            
            %define the normalized transition rate matrix
            normedrates = transit_rates./sigma_c_km;
            normedrates_vec = normedrates(:);
            
            % find the smallest j satisfies the (sigma ckm from 1 to j) <=
            % U2 * (sigma ckm from m =1 and m~=k)
            i_select=min(find(cumsum(normedrates_vec)>=U2));
            
            area_mat_vec = area_mat(:); % stack the matrix to vector
            
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
            area_mat = reshape(area_mat_vec,L,L); % reshape vector into matrix back
        
        % update t_exit and store the stats at this iteration
        t_exit = t_exit  + tau;
        N_S(k,  t_counter) = NS;
        N_I(k,  t_counter) = NI;
        N_R(k,  t_counter) = NR;
        t_counter = t_counter +1;
        
        %%%%%%%%%%%%% uncomment the code here so that you can watch a moive : )
        
        %set(h11,'CData',area_mat)
        %keyboard
        %pause(0.01);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        end
  
t_exit_mat(k,1)=t_exit;
k = k+1;
figure;
imagesc(area_mat)
mytitle1=['Snapshot when r remove =  ',num2str(r_remove)];
title(mytitle1);
end
toc