%%%%% Homework 3 (c)%%%%%%%%%
%%%%% author: Mingfu Liang %%%%%%%%
%%%%% date: 03/05/2019 %%%%%%%

tic
p_immune_count = 1;
range = 0.4/0.02 +1;
mean_t_exit = zeros(range,1);
std_t_exit = zeros(range,1);
mean_ND = zeros(range,1);
std_ND = zeros(range,1);

for p_immune = 0%:0.02:0.4

text = ['p_immune = ', num2str(p_immune), ' \n'];
fprintf(text);
    
L =100;
L_2 = L*L;
N_trials = 100;
r_inject = 0.5;
r_remove = 1.1;
N_S = zeros(N_trials,1);
N_I = zeros(N_trials,1);
N_R = zeros(N_trials,1);
t_exit_mat = zeros(N_trials,1);
N_D_mat =zeros(N_trials,1);
posit = 1: L; % define the index variables 
up_shift = circshift(posit,1); % shift the variables up one unit
down_shift = circshift(posit,-1); % shift the variables down one unit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i =1: N_trials
    
    % Starting with an initial condition in which only a single individual in the center of the
    % system (i = 10; j = 10) is infected (L = 20)
    % here denote 1 is susceptible, 2 is infected, 3 is removed
    
    immune_num = L_2 * p_immune;
    area_mat = ones(L, L);
    
    for select_ind = 1:immune_num
              i_select = i_ind(select_ind);
              j_select = j_ind(select_ind);
              area_mat(i_select, j_select) =3;
    end
    
    area_mat(50,50) = 2;
    t_counter =1;
    t_exit = 0;
    NS = L_2 - immune_num;
    NI = 1;
    NR = immune_num;
    
    %transit_rates = zeros(L, L);
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
            U1 = rand();
            tau = -log(U1)/sigma_c_km;
            U2 = rand();
            sigma_ckm = 0;
            normedrates = transit_rates./sigma_c_km;
            normedrates_vec = normedrates(:);
            i_select=find(cumsum(normedrates_vec)>=U2, 1 );
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
    
        N_S(i,  t_counter) = NS;
        N_I(i,  t_counter) = NI;
        N_R(i,  t_counter) = NR;
        t_counter = t_counter +1;
        end
        
    t_exit_mat(i,1)=t_exit;
    N_D = N_R(i,t_counter-1) - immune_num;
    N_D_mat(i,1) = N_D;
    
    %set(h11,'CData',area_mat)
    %pause(0.01);
end

mean_t_exit(p_immune_count,1)=mean(t_exit_mat);
std_t_exit(p_immune_count,1)=std(t_exit_mat);
mean_ND(p_immune_count,1)=mean(N_D_mat);
std_ND(p_immune_count,1)=std(N_D_mat);
p_immune_count = p_immune_count + 1;
end
toc

t_range = 0:0.02:0.4;
figure;
plot(t_range,mean_t_exit);
title2=['mean of t exit verse p immue'];
title(title2)
xlabel('$p_{immune}$','Interpreter','latex','FontSize',13)
ylabel('$<t_{exit}>$','Interpreter','latex','FontSize',13)

figure;
plot(t_range,std_t_exit);
xlabel('$p_{immune}$','Interpreter','latex','FontSize',13)
ylabel('$std(t_{exit})$','Interpreter','latex','FontSize',13)
title3=['standard deviation of t exit verse p immue'];
title(title3)

figure;
plot(t_range,mean_ND);
xlabel('$p_{immune}$','Interpreter','latex','FontSize',13)
ylabel('$<N_{D}>$','Interpreter','latex','FontSize',13)
title4=['mean of ND verse p immue'];
title(title4)

figure;
plot(t_range,std_ND);
xlabel('$p_{immune}$','Interpreter','latex','FontSize',13)
ylabel('$std(N_{D})$','Interpreter','latex','FontSize',13)
title5=['standard deviation of ND verse p immue'];
title(title5)