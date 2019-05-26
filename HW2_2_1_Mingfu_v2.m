%%%%%%% author: Mingfu Liang
%%%%%%% date: 02/27/2019
%%%%%%% Problem 2 Phase Transition in the Ising Model
tic

%%%%%%%%%% initialize the parameter %%%%%%%%%%%%%%%%%%%
H=0; % H is non-dimensionalized version of B
J=1; % J is non-dimensionalized coupling strength
L=25; % Consider a two-dimensional system of spin 1/2 magnetic dipoles arrayed on a square L*L
L_2 = L^2; % The size of the one Monte-Carlo step, L_2 is the number of spins in the system
t_corr = 200; % fisrt establish equilibrium by taking t_corr Monte-Carlo steps as intial configuration
N_trials = 400; % take N_trials * t_corr Monte-Carlo steps
T_high = 3; % The highest value in the relevant range of the temperature
T_low =1.5; % The lowest value in the relevant range of the temperature
delta_T = -0.1; % The change of temperature in each steps
T_range =zeros(16,1);% Initialize the range of temperature


% Generate the time range which should be [1.5,3.0]
h=1;
for T = T_high:delta_T:T_low
    T_range(h,1)=T;
    h=h+1;
end

posit = 1: L; % define the index variables 
up_shift = circshift(posit,1); % shift the variables up one unit
down_shift = circshift(posit,-1); % shift the variables down one unit
counter =0; % count the total iteration, using for debug 
T_init = T_high; % To calculate the initialization configuration

%%%%%%%%%% Initialize the M_{Sigma}, U_{Sigma} and the flip counter matrix
M_Sigma= zeros(16,400);
U_Sigma_matrix=zeros(16,400);
U_Sigma = zeros(16,1);
U_Sigma_Square =zeros(16,400);
flip_tol = zeros(16,1);


%%%%%%%%%% Metropolis Algorithm using the Ising model %%%%%%%

%%%%%%%%%% Then take N_trials*t_corr Monte-Carlo steps and measure
%%%%%%%%%% the expectation values listed below by sampling the corresponding quantities every
%%%%%%%%%% t_corr Monte-Carlo steps.
temp_count = 1;
for T = T_high:delta_T:T_low
    %%%%%%%%%% first establish equilibrium by taking t_corr Monte-Carlo
    %%%%%%%%%% steps using an initial configuration
    
            if exist('Sigma_mat')==0
                Sigma_mat = randsrc(L); % generate square matrix L*L and each lattice point (i,j) the dipole is only 1 or -1
            end
            for k = 1:1:t_corr
                    [i_ind,j_ind]=ind2sub([L,L],randperm(L^2)); % in each Monte-Carlo step generate a random permutation of all spin indices
                for select_ind = 1: L_2 % perform one Monte-Carto step with L*L attempted spin flips
                i_select = i_ind(select_ind);
                j_select = j_ind(select_ind);
                Sigma_select_old = Sigma_mat(i_select, j_select);
                Sigma_select_new = -1*Sigma_mat(i_select, j_select);% propose the candidate spin flip
                Delta_U = 2*H*Sigma_select_old + 2*J*Sigma_select_old*(Sigma_mat(up_shift(i_select), j_select)+Sigma_mat(down_shift(i_select), j_select)+Sigma_mat(i_select, up_shift(j_select))+Sigma_mat(i_select, down_shift(j_select)));
                accept_prob = min(exp(-Delta_U/T),1);
                U = rand;
                        if U < accept_prob
                            Sigma_mat(i_select,j_select)=Sigma_select_new;
                        end
                end
            end    
    
        % initialize the accepted flip counter
        flip_count =0;
        
        for N=1:N_trials
            
            for t = 1:t_corr
                [i_ind,j_ind]=ind2sub([L,L],randperm(L^2)); % in each Monte-Carlo step generate a random permutation of all spin indices
                for select_ind = 1: L_2 % perform one Monte-Carto step with L*L attempted spin flips
                        i_select = i_ind(select_ind);
                        j_select = j_ind(select_ind);
                        % define the Sigma_{select_old} for energy
                        % difference calculation
                        Sigma_select_old = Sigma_mat(i_select,j_select);
                        % propose the candidate spin flip
                        Sigma_select_new = -1*Sigma_mat(i_select,j_select);
                        % Compute the energy difference Delta_{U} using (5.3)
                        Delta_U = 2*H*Sigma_select_old+2*J*Sigma_select_old*(Sigma_mat(up_shift(i_select),j_select)+Sigma_mat(down_shift(i_select),j_select)+Sigma_mat(i_select,up_shift(j_select))+Sigma_mat(i_select,down_shift(j_select)));
                        accept_prob = min(exp(-Delta_U/T),1);
                        U = rand;
                        if U < accept_prob
                            Sigma_mat(i_select,j_select)=Sigma_select_new;
                            flip_count = flip_count+1;
                        end
                end
                counter = counter+1;
            end
            
            % Calculate M_Sigma and U_Sigma
            for select_ind = 1: L_2 
            i_select = i_ind(select_ind);
            j_select = j_ind(select_ind);
            M_Sigma(temp_count,N)=M_Sigma(temp_count,N)+ Sigma_mat(i_select,j_select);
            U_Sigma(temp_count)=U_Sigma(temp_count)-J*Sigma_mat(i_select,j_select)*(Sigma_mat(up_shift(i_select),j_select)+Sigma_mat(down_shift(i_select),j_select)+Sigma_mat(i_select,up_shift(j_select))+Sigma_mat(i_select,down_shift(j_select)))/2;
            U_Sigma_Square(temp_count,N)=U_Sigma_Square(temp_count,N)+(-J*Sigma_mat(i_select,j_select)*(Sigma_mat(up_shift(i_select),j_select)+Sigma_mat(down_shift(i_select),j_select)+Sigma_mat(i_select,up_shift(j_select))+Sigma_mat(i_select,down_shift(j_select)))/2);
            %U_Sigma_matrix(temp_count,N)=-J*Sigma_mat(i_select,j_select)*(Sigma_mat(up_shift(i_select),j_select)+Sigma_mat(down_shift(i_select),j_select)+Sigma_mat(i_select,up_shift(j_select))+Sigma_mat(i_select,down_shift(j_select)))/2;
            end
            
        end
      flip_tol(temp_count,1)= flip_count;
      temp_count=temp_count+1; 
      
      %%%%%%%%%% Snapchat
        if T == 2
            figure;
            imagesc(Sigma_mat)
            colormap(gray(2))
            axis ij
            axis square
            mytitle5=['Temperature plot when Temperature = 2.0, L = ',num2str(L)];
            title(mytitle5);
        end
        if T == 2.5
            figure;
            imagesc(Sigma_mat)
            colormap(gray(2))
            axis ij
            axis square
            mytitle6=['Temperature plot when Temperature = 2.5, L = ',num2str(L)];
            title(mytitle6)
        end
        if T == 3
            figure;
            imagesc(Sigma_mat)
            colormap(gray(2))
            axis ij
            axis square
            mytitle11=['Temperature plot when Temperature = 3, L = ',num2str(L)];
            title(mytitle11)
        end
end

%%%%%%%%%% Equilibrium Temperature plot when Temperature = 1.5 %%%%%%%%%
figure;
imagesc(Sigma_mat)
colormap(gray(2))
axis ij
axis square
mytitle8=['Equilibrium Temperature plot when Temperature = 1.5, L = ',num2str(L)];
title(mytitle8)

%%%%%%%%%% measure the expectation values listed below by sampling the
%%%%%%%%%% corresponding quantities 

%%% Fraction of accepted spin flips, here flip_tol() is the total number of
%%% spin flips of different temperature, which is 16*1 vector, then divide
%%% by the total attempt of spin flips
Filp_Faction_Accp = flip_tol()/(L_2*t_corr*N_trials);

%%% Define the M(Sigma), which is 16*400 matrix, each row has 400 sample
%%% obtained from specified temperature
%M_Sigma_vec = M_Sigma/L_2;

%%% Define U(T_tiuda) from samples
U_T_tiuda = U_Sigma/(L_2*N_trials);

%%% Define the Magnetization of whole state from samples
M_T_tiuda =sum(abs(M_Sigma),2)/N_trials;

%%% Define the Magnetization per spin from samples
m_T_tiuda = sum(abs(M_Sigma),2)/(N_trials*L_2);

%%% Define the <U(\tiuda(T))^2> from samples
U_T_tiuda_square= sum(U_Sigma_Square.^2,2)/(N_trials);

%%% Define the <M(\Sigma)^2> from samples
M_T_tiuda_square = sum(M_Sigma.^2,2)/N_trials;

%%% Define Sepcific Heat from samples
SepcificHeat = (U_T_tiuda_square - ((U_Sigma/N_trials).^2))./(L_2*T_range.^2);

%%% Define Susceptibility without divide L^2 from samples
Susceptibility = (M_T_tiuda_square-M_T_tiuda.^2)./T_range;

%%% Define the correct Susceptibility from samples
Susceptibility_L_2 = (M_T_tiuda_square-M_T_tiuda.^2)./(T_range*L_2);

%%% Visualize the statistics
figure;
plot(T_range,Susceptibility_L_2,'linewidth',2)
mytitle14 =['Susceptibility with the new equation which divide L^2 when L is ',num2str(L)];
title(mytitle14);
figure;
plot(T_range,U_T_tiuda,'linewidth',2)
mytitle1 =['Energy per spin when L is ',num2str(L)];
title(mytitle1);
figure;
plot(T_range,SepcificHeat,'linewidth',2)
mytitle2 =['Specific Heat when L is ',num2str(L)];
title(mytitle2);
figure;
plot(T_range,Susceptibility,'linewidth',2)
mytitle3 =['Susceptibility of whole state when L is ',num2str(L)];
title(mytitle3);
figure;
plot(T_range,m_T_tiuda,'linewidth',2)
mytitle4 =['Magnetization per spin when L is ',num2str(L)];
title(mytitle4);
figure;
plot(T_range,Filp_Faction_Accp,'linewidth',2)
mytitle10 =['Fraction of accepted spin flips when L is ',num2str(L)];
title(mytitle10);
toc