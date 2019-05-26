%%%%% Homework 3 a %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/05/2019 %%%%%%%
tic

L =20;
L_2 = L*L;
N_trials = 1000;
r_inject = 0.5;
r_remove = 1.1;
N_S = zeros(1000,1);
N_I = zeros(1000,1);
N_R = zeros(1000,1);
t_exit_mat = zeros(1000,1);
posit = 1: L; % define the index variables 
up_shift = circshift(posit,1); % shift the variables up one unit
down_shift = circshift(posit,-1); % shift the variables down one unit

%%%%%%%%%%%

for i =1: N_trials
    
    % Starting with an initial condition in which only a single individual in the center of the
    % system (i = 10; j = 10) is infected (L = 20)
    % here denote 1 is susceptible, 2 is infected, 3 is removed
    
    area_mat = ones(L, L);
    area_mat(10,10) = 2;
    t_counter =1;
    t_exit = 0;
    transit_rates = zeros(L, L);
    ckm_mat = zeros(L, L);
    [i_ind,j_ind]=ind2sub([L, L],randperm(L^2)); % generate a random permutation of all index
    %figure;
    %h11=imagesc(area_mat);
    NS = 399;
    NI = 1;
    NR = 0;
    
        while NI~=0
            for select_ind = 1: L_2
        
                    i_select = i_ind(select_ind);
                    j_select = j_ind(select_ind);
                if area_mat(i_select, j_select) ==1 % which means that this lattice point is  susceptible
              
                transit_rates(i_select, j_select) = r_inject*( 0 + (area_mat(up_shift(i_select), j_select)==2) ...
                                                                          +(area_mat(down_shift(i_select), j_select)==2) ...
                                                                             +(area_mat(i_select, up_shift(j_select)) ==2)...
                                                                          +(area_mat(i_select, down_shift(j_select)) ==2) ...
                                                            +(area_mat(up_shift(i_select), down_shift(j_select)) ==2) ...
                                                         +(area_mat(down_shift(i_select), down_shift(j_select))==2) ...
                                                                   +(area_mat(up_shift(i_select), up_shift(j_select)) ==2) ...
                                                             +(area_mat(down_shift(i_select), up_shift(j_select)) ==2));
                end
          
                if area_mat(i_select, j_select) == 2 % which means that this lattice point is infected
                        transit_rates(i_select, j_select) = r_remove;
                end
                
            end

            
            sigma_c_km = sum(sum(transit_rates,2));
    
            U1 = rand();
            tau = -log(U1)/sigma_c_km;
            U2 = rand();
            sigma_ckm = 0;
    
        for select_ind = 1: L_2
        
          i_select = i_ind(select_ind);
          j_select = j_ind(select_ind);
          
          sigma_ckm = sigma_ckm + transit_rates(i_select, j_select);
          
            if sigma_ckm >=U2*sigma_c_km
            %%% change the state
            t_exit = t_exit  + tau;
                if area_mat(i_select, j_select) ~=3
                    area_mat(i_select, j_select) = area_mat(i_select, j_select) + 1;
                end
                break
            end
        end
    
        S = find(area_mat == 1);
        I = find(area_mat == 2);
        R = find(area_mat == 3);
    
        NS = size(S,1);
        NI  = size(I,1);
        NR = size(R,1);
    
    
        N_S(i,  t_counter) = NS;
        N_I(i,  t_counter) = NI;
        N_R(i,  t_counter) = NR;
        t_counter = t_counter +1;
        end
        
    t_exit_mat(i,1)=t_exit;
    %set(h11,'CData',area_mat)
    pause(0.01);
end

hist(t_exit_mat)
mean(t_exit_mat)
toc
