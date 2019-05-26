%%%%% Homework 3%%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/05/2019 %%%%%%%
function MingfuLiang_HW3(n_trial,L)
tic


L_2 = L*L;
N_trials = n_trial;
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
            Infected_Neigh=area_mat==2;
             %%% calculate I_neighbor for sus
            [i_col,j_col] = find(area_mat == 1);
            for select_ind = 1:size(i_col,1)
        
              i_select = i_col(select_ind);
              j_select = j_col(select_ind);
              
              transit_rates(i_select, j_select) = r_inject*( 0 + Infected_Neigh(up_shift(i_select), j_select) ...
                                                                          +Infected_Neigh(down_shift(i_select), j_select) ...
                                                                             +Infected_Neigh(i_select, up_shift(j_select)) ...
                                                                          +Infected_Neigh(i_select, down_shift(j_select))  ...
                                                            +Infected_Neigh(up_shift(i_select), down_shift(j_select)) ...
                                                         +Infected_Neigh(down_shift(i_select), down_shift(j_select)) ...
                                                                   +Infected_Neigh(up_shift(i_select), up_shift(j_select))  ...
                                                             +Infected_Neigh(down_shift(i_select), up_shift(j_select)));


                
            end
            transit_rates(area_mat==2) =r_remove;
            sigma_c_km= sum(transit_rates(:));
    
            U1 = rand();
            tau = -log(U1)/sigma_c_km;
            U2 = rand();
            sigma_ckm = 0;
            
%             normedrates = transit_rates./sigma_c_km;
%             normedrates_vec = normedrates(:);
%             i_select=min(find(cumsum(normedrates_vec)>=U2));
%             area_mat_vec = area_mat(:);
%             if area_mat_vec(i_select) ==2
%                 area_mat_vec(i_select) =3;
%                 NR = NR +1;
%                 NI = NI - 1;
%             elseif area_mat_vec(i_select) ==1
%                 area_mat_vec(i_select) = 2;
%                 NI = NI + 1;
%                 NS = NS - 1;
%             end
%             area_mat = reshape(area_mat_vec,L,L)';
            
        for select_ind = 1: L_2
        
          i_select = i_ind(select_ind);
          j_select = j_ind(select_ind);
          
          sigma_ckm = sigma_ckm + transit_rates(i_select, j_select);
          
            if sigma_ckm >=U2*sigma_c_km
                if area_mat(i_select, j_select)==1
                     area_mat(i_select, j_select) = 2;
                     NI = NI + 1;
                     NS = NS - 1;
                 elseif area_mat(i_select, j_select)==2
                      area_mat(i_select, j_select) = 3;
                      NR = NR +1;
                      NI = NI - 1;
                 end
                break
            end
        end
                
               
            t_exit = t_exit  + tau;
            
%         S = find(area_mat == 1);
%         I = find(area_mat == 2);
%         R = find(area_mat == 3);
%     
%         NS = size(S,1);
%         NI  = size(I,1);
%         NR = size(R,1);
    
    
        N_S(i,  t_counter) = NS;
        N_I(i,  t_counter) = NI;
        N_R(i,  t_counter) = NR;
        t_counter = t_counter +1;
        end
        
    t_exit_mat(i,1)=t_exit;
    %set(h11,'CData',area_mat)
    %pause(0.01);
end

hist(t_exit_mat)
mean(t_exit_mat)
toc
end
