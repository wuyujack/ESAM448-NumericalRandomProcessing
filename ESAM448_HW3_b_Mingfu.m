%%%%% Homework 3 b %%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/05/2019 %%%%%%%

tic
L =100;
L_2 = L*L;
N_trials = 1;
r_inject = 0.5;
N_S = zeros(6,1);
N_I = zeros(6,1);
N_R = zeros(6,1);
t_exit_mat = zeros(6,1);
[i_ind,j_ind]=ind2sub([L, L],randperm(L^2)); % generate a random permutation of all index
k=1;
for r_remove = 0:0.3:1.5
    

%r_remove = 1.1;

posit = 1: L; % define the index variables 
up_shift = circshift(posit,1); % shift the variables up one unit
down_shift = circshift(posit,-1); % shift the variables down one unit

%%%%%%%%%%%


    
    % Starting with an initial condition in which only a single individual in the center of the
    % system (i = 10; j = 10) is infected (L = 20)
    % here denote 1 is susceptible, 2 is infected, 3 is removed
    
    area_mat = ones(L, L);
    area_mat(50,47:53) = 2;
    t_counter =1;
    t_exit = 0;
    transit_rates = zeros(L, L);
    
%     figure;
%     h11=imagesc(area_mat);
%     title1 = ['snapshot of r remove = ',num2str(r_remove)];
%     title(title1)
    
    NS = 10000-7;
    NI = 7;
    NR = 0;
    
        while NI~=0
            if NI == 10000
                break
            end
            Infected_Neigh=area_mat==2;
            for select_ind = 1: L_2
        
              i_select = i_ind(select_ind);
              j_select = j_ind(select_ind);
              
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
                    %area_mat(i_select, j_select) = area_mat(i_select, j_select) + 1;
                end
                break
            end
        end
                 if area_mat(i_select, j_select)==1
                     area_mat(i_select, j_select) = 2;
                     NI = NI + 1;
                     NS = NS - 1;
                 elseif area_mat(i_select, j_select)==2
                      area_mat(i_select, j_select) = 3;
                      NR = NR +1;
                      NI = NI - 1;
                 end
%         S = find(area_mat == 1);
%         I = find(area_mat == 2);
%         R = find(area_mat == 3);
%     
%         NS = size(S,1);
%         NI  = size(I,1);
%         NR = size(R,1);
    
    
        N_S(k,  t_counter) = NS;
        N_I(k,  t_counter) = NI;
        N_R(k,  t_counter) = NR;
        t_counter = t_counter +1;
        %set(h11,'CData',area_mat)
        %pause(0.0000001);
        end
  
t_exit_mat(k,1)=t_exit;
k = k+1;
%figure;
%imagesc(area_mat)
%mytitle1=['Snapshot when r remove =  ',num2str(r_remove)];
%title(mytitle1);
end

%figure;
%hist(t_exit_mat)
%mean(t_exit_mat)
figure;
r_remove = 0:0.3:1.5;
plot(r_remove,t_exit_mat);
toc