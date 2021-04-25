% ODE model for Drosophila suzukii
%Initialization
clear, clc % To clear the workspace and the command function

%% PARAMETERS functions

%Growth rate (Drosophila suzukii)
 a = 1.2*(10^(-4));
 T_l = 3;
 T_m = 30;
 m = 6;

%Mortality rate (Drosophila suzukii)
a1 = -5.4E-06;
b1 = 0.0005194;
c1 = -0.0116827;
d1 = 2.16E-05;
e1 = 1.3146586;

% Birth rate (Drosophila suzukii)
alpha = 659.06;
gamma = 88.53;
lambda = 52.32;
delta = 6.06;
tau = 22.87;

% Sex ratio
S_R = 0.5;

% Mating ratio
R_remate = 0;
R_mate = 1;

% Additional parameters
N_stages = 8; % Eggs/3 larva stages/ pupa /male/ unmated female/ mated female

%% DATA

%Load additional data. To update depending on the temperature used
Temp_avg = csvread('TemperaturesF1.csv', 1);

%Compute the rates based on temperature
G_R = growth_rate(Temp_avg(1),a,T_l,T_m,m); %Calling the growth rate function

B_R = birth_rate_suzuki(Temp_avg(1),alpha,gamma,lambda,delta,tau);%Calling the birth rate function

M_R = mortality_rate(Temp_avg(1),a1,b1,c1,d1,e1); %Calling the mortality rate function

%Initialize stages
Pest_stages(1:N_stages) = stage; % We create an array of the stage class
Pest_stages = Initialize_stages_ode(B_R,M_R,G_R,S_R,R_mate,R_remate,Pest_stages); %We initialize the parameters associated to each stage class


%% DYNAMIC MODEL

%Initial values for the variables
x = zeros(N_stages,1); %Array of variables

A_cont =compute_A_continous(Pest_stages); %We build the continuous time A matrix
% Discretization of the system based on the zero order holder
sysc = ss(A_cont,[],eye(8),[]);
sysd = c2d(sysc,1,'zoh');
A_dis = sysd.A; %Discretized matrix A

%% SIMULATIONS
Simulation_time =length(Temp_avg); %Simulation lenght based on the temperature array introduced

%Initial conditions
x(1) = 10000000; % eggs
x(7) = 0; % adult males
x(8) = 10000000; % adult mated females

x_hist=x; %We start the array x_hist to keep a historic of the evolution of the states

for t=1:(Simulation_time-1) %Loop for the simulation
    
    x = A_dis*x; %Compute the state at the next time step
 
    x_hist = [x_hist,x]; %Store the states
    
       
    % First recompute rates
    G_R = growth_rate(Temp_avg(t+1),a,T_l,T_m,m);
    B_R = birth_rate_suzuki(Temp_avg(1),alpha,gamma,lambda,delta,tau);
    M_R = mortality_rate(Temp_avg(t+1),a1,b1,c1,d1,e1);
    Pest_stages = Initialize_stages_ode(B_R,M_R,G_R,S_R,R_mate,R_remate,Pest_stages);
   
    %Update the matrix A
    A_cont =compute_A_continous(Pest_stages);
    sysc = ss(A_cont,[],eye(8),[]);
    sysd = c2d(sysc,1,'zoh');
    A_dis = sysd.A;
    
end


%% PLOTS
t1 = datetime(2017,3,15,12,0,0);
t=t1+days(0:Simulation_time-1);

tt = 1:Simulation_time;

% Read the validation data

valid_1 = csvread('ValidField1.csv', 1);

% Plot adult males population

figure
plot(t,x_hist(6,:),'LineWidth',2);
hold on
scatter(t(valid_1(:,1)),valid_1(:,2),'LineWidth',2)
hold on
plot(t(valid_1(:,1)),valid_1(:,2),":",'LineWidth',1.5)
legend("Adult males", "Field data",'Fontsize',15);
xlabel('Time [days]','Fontsize',15);
ylabel('Number of adult males','Fontsize',15);
title('Montelibretti Field 2017','Fontsize',20);
set(gca,'FontSize',15)

% Write data on a txt file

Tab = table(transpose(tt), transpose(x_hist));
writetable(Tab, 'Allstages.csv', 'Delimiter', ';');
