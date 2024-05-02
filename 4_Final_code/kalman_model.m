function [n_states_signal, phi_kalman_model, C_kalman_model, B_kalman_model, tau_s, C_signal]=kalman_model(tau_tank)

set(0, 'DefaultLineLineWidth',2);
set(0, 'DefaultaxesLineWidth',1);
set(0, 'DefaultaxesFontSize',12);
set(0, 'DefaultTextFontSize',12);
set(0, 'DefaultAxesFontName','Times');


%% Generate demand curve

generate_demad_curve
freq_sampling=1;

%% Finding model for the demand curve

n_demand_points=24*60;  % Number of data points for Fourier analysis

%   Amplitude spectrum plot

Y = fft(D_demand,n_demand_points);
P2 = abs(Y/n_demand_points);
P1 = P2(1:n_demand_points/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = freq_sampling*(0:(n_demand_points/2))/n_demand_points;


figure
plot(t,-D_demand,'b')
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend('Consumer demand curve (-ve)');
title('Total consumer demand');
xlim([0 T_min]);
xticks(0:2:T_min)
grid

figure
stem(f,P1,'b','LineWidth',2)
xlabel('Frequency [Hz]');ylabel('|P(f)|')
title('Amplitude spectrum of the consumer demand curve');
xlim([0 0.005])
xticks(0:0.0005:0.005)
grid

% Determine dominant frequencies

f_cutoff=0.02;   % Cut-off amplitude

f_index=find(P1>f_cutoff);     % Selecting frequencies above cut-off amplitude
f_set=nonzeros(f(f_index));
n_f_set=length(f_set);

% Finding parameters of the curve using least square

f_omega=2*pi*f_set;

n_D_demand=length(D_demand);    % number of data points in demand data

n=(1:n_D_demand)';

A0=ones(n_D_demand,1);          % A1 matrix for a_0 coeff
A=cos(n*f_omega');               % A2 matrix for a_n coeff corresponding to cos() terms
B_kalman_model=sin(n*f_omega');               % A3 matrix for b_n coeff corresponding to sin() terms

M=[A0 A B_kalman_model];                   % A matrix

coeff_curve_model=(M'*M)\M'*D_demand';      % Coeff of the model in order [a_0;a_1;...;a_n;b_1;...;b_n]

a_0=coeff_curve_model(1);
a=coeff_curve_model(2:1+n_f_set);
b=coeff_curve_model(2+n_f_set:1+2*n_f_set);

%% Recreating the signal using state space model (Open Loop)

n_states_signal=2*n_f_set+1;   % number of states equal to 2*number of frequency + a_0

tau_s=60;                % Sampling time

% Creating the phi matrix for signal recreation SS model
phi_signal=zeros(n_states_signal);
C_signal=zeros(1,n_states_signal);

phi_signal(1,1)=1;
C_signal(1,1)=1;

for i=1:n_f_set
    phi_signal(2*i:2*i+1,2*i:2*i+1)=[cos(f_omega(i)*tau_s) sin(f_omega(i)*tau_s);...
        -sin(f_omega(i)*tau_s) cos(f_omega(i)*tau_s)];
    
    C_signal(1,2*i:2*i+1)=[1 0];
end

Obsv_signal=obsv(phi_signal,C_signal);
Obsv_signal_cond=cond(Obsv_signal);

T_sim=T_min*60/tau_s;     % Simulation time for signal recreation

% Initialize the states at k=1
x_rec=a_0;
for i=1:n_f_set
    x_rec=[x_rec;a(i)*cos(0)+b(i)*sin(0); -a(i)*sin(0)+b(i)*cos(0)];
end

x_rec=[x_rec,zeros(n_states_signal,T_sim)];
y_rec=zeros(1,T_sim);

% Running simulation for signal recreation

for k=1:T_sim
    
    x_rec(:,k+1)=phi_signal*x_rec(:,k);
    y_rec(:,k)=C_signal*x_rec(:,k);
    
end
y_rec(:,k+1)=C_signal*x_rec(:,k+1);
D_recreate=y_rec;

time_plot=linspace(0,T_min,length(D_recreate));

figure
plot(t,-D_demand,'b',time_plot,-D_recreate,'g:')
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend('Actual consumer demand curve (-ve)','Predicted consumer demand curve (-ve)');
title('Comparison of actual and predicted consumer demand curve from state space model');
xlim([0 T_min]);
xticks(0:2:T_min)
grid


%% Predicition with Kalman filter

% Creating state space model for Kalman filter
phi_kalman_model=[phi_signal,zeros(n_states_signal,1)];
phi_kalman_model=[phi_kalman_model;[tau_tank*tau_s,tau_tank*tau_s,0,tau_tank*tau_s,0,1]];    % phi matrix with augmented state
phi_kalman_model(end,end)=1;

C_kalman_model=[zeros(1,n_states_signal) 1];  % Output matrix C

Obsv_kalman_model=obsv(phi_kalman_model,C_kalman_model);  % observability matrix
cond_Obsv=cond(Obsv_kalman_model);           % condtion number of observability matrix
rank_Obs_matrix=rank(Obsv_kalman_model);     % rank of observability matrix
n_states_signal+1;                             % number of states

B_kalman_model=[zeros(n_states_signal,1);tau_tank*tau_s];
