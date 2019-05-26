tic

L = 20;
n_trial = 1000;
L_2 = L*L;

N_trials = n_trial;
r_inject = 0.5;
r_remove = 1.1;
N_S = zeros(1000,400);
N_I = zeros(1000,400);
N_R = zeros(1000,400);
t_exit_mat = zeros(1000,1);
posit = 1: L; % define the index variables 
up_shift = circshift(posit,1); % shift the variables up one unit
down_shift = circshift(posit,-1); % shift the variables down one unit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i =1: N_trials
    
    % Starting with an initial condition in which only a single individual in the center of the
    % system (i = 10; j = 10) is infected (L = 20)
    % here denote 1 is susceptible, 2 is infected, 3 is removed
    
    area_mat = ones(L, L);
    area_mat(10,10) = 2;
    t_counter =1;
    t_exit = 0;
    NS = 399;
    NI = 1;
    NR = 0;
    
    %transit_rates = zeros(L, L);
    %ckm_mat = zeros(L, L);
    %figure;
    %h11=imagesc(area_mat);

        while NI~=0
                       Infected_Neigh=area_mat==2;
                       
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

            sigma_c_km= sum(transit_rates(:));
    
            
            tau = -log(rand())/sigma_c_km;
            
            normedrates = transit_rates/sigma_c_km;
            normedrates_vec = normedrates(:);
            i_select=find(cumsum(normedrates_vec)>=rand(), 1);
            area_mat_vec = area_mat(:);
            if area_mat_vec(i_select) ==2
                area_mat_vec(i_select) =3;
                NR = NR +1;
                NI = NI - 1;
            elseif area_mat_vec(i_select) ==1
                area_mat_vec(i_select) = 2;
                NI = NI + 1;
                NS = NS - 1;
            end
           area_mat = reshape(area_mat_vec,L,L);
               
            t_exit = t_exit  + tau;
            
        N_S(i,  t_counter) = NS;
        N_I(i,  t_counter) = NI;
        N_R(i,  t_counter) = NR;
        t_counter = t_counter +1;
        end
        
    t_exit_mat(i,1)=t_exit;
    %set(h11,'CData',area_mat)
    %pause(0.01);
end

hist(t_exit_mat,20);
mean_t_exit = mean(t_exit_mat);
text = ['mean of t_exit is ', num2str(mean_t_exit),' \n'];
fprintf(text);
t_exit = t_exit_mat;
toc