%% controller_predictor


%% Simulation Parameters

tic

ts=1;           % [s] Sampling time for the model
Ts=1;             % [min] Sampling time for MPC
T_min=24*6;         % [min] Simulation time. Provide simulation time in multiple of 24 only.
Hp=24;            % [min instance] Prediction horizon

% Simulation

T_sec=T_min*60;      % [s] Simulation time in seconds

N=(T_sec+(Hp*60))/ts;         % [] Simulation instances

% Generating demand curve for simulation. Number of data points to be for Hp(24 min) more than the
% simulation time. Data generated with sampling time of 1 sec

t=linspace(0,T_min+Hp,(N+1));
d2=(-0.3*(3*sin((t/(3.83))+5)-4*cos(2*t/(3.83))+8)*0.04)-0.02; % [m^3/s] -ve-> flow out of network
d5=(-0.18*(3*sin((t/(3.83))+5)-4*cos(2*t/(3.83))+8)*0.04)-0.02; % [m^3/s] -ve-> flow out of network

v1=mean(d2./(d2+d5));
v2=mean(d5./(d2+d5));

std_process_noise=0.005;

d2=d2+randn(1,(N+1))*std_process_noise;
d5=d5+randn(1,(N+1))*std_process_noise;

dc=[d2;d5];  % [m^3/h] Flow demand at consumer node


% Generating cost of electricity

C_con_set=[0.7*ones(1,6),1.4*ones(1,14),0.7*ones(1,4)];

% C_con_set=[1*ones(1,6),1*ones(1,14),1*ones(1,4)];

% C_con_set=[0.7,0.7,0.7,0.7,0.8,0.8...
%     0.9,1,1.2,1.4,1.35,1.3,1.25,1.2,1.25,1.3,1.35,1.4,1.2,1.1,...
%     0.9,0.7,0.7,0.7];

C_con=[];
for i=1:(T_min+Hp)/24
    C_con=horzcat(C_con,C_con_set);
end

p_0_all=zeros(1,(Ts*60/ts)+1);
prompt = 'Please provide initial level in the tank [m] ';
p_0_all(:,1)= input(prompt)*rho_fluid*g/1e5;   % Ask for the reference node

std_meas_noise=0.0005;

disp(' ');
U_con=[];
Q_C=[];
Q_T=[];
P_BAR=[];
D_TAU=[];
P_val_0_noise=p_0_all(:,1)+randn*std_meas_noise;
EPS_VAL=[];
P_val_0=p_0_all(:,1);
P_val_esti=[];

NMPC_problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
Pp_kalman=100*diag(1*ones(n_states_signal+1,1));               %Autocovariance for the initial estimate
Q_kalman=(std_process_noise^2)*2*diag(ones(n_states_signal+1,1));%std_state_noise^2;                %Covariance of system noise
Q_kalman(end,end)=Q_kalman(end,end)*(tau*tau_s)^2;
R_kalman=std_meas_noise^2;                    %Covariance of measurement noise

x_hat=zeros(6,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j=1;

for i=1:T_min
    i
    
    % Extracting cost and demand for next Hp instances
    C_Hp=C_con(:,i:i+Hp-1);
    %     dc_Hp=dc_con(:,i:i+Hp-1);

    p_0_noise=P_val_0_noise(:,end);
    
    predictor_NMPC

    D_predict(D_predict>0)=0;
    D_predict(D_predict<-0.45)=-0.45;

    dc_Hp=[v1*D_predict;v2*D_predict];

    P_val_esti=horzcat(P_val_esti,p_esti);


    opti.set_value(p_0_par,p_esti);
    opti.set_value(dc_par,dc_Hp);
    opti.set_value(C_par,C_Hp);
    
    
    sol = opti.solve();
    u_all= sol.value(u_var);
    eps_all=sol.value(eps_var);
    
    u_con=u_all(:,1);
    eps_con=eps_all(:,1);
    
    
    %Predict state for next time step
    x_hat(:,i+1)= phi_kalman_model*x_hat(:,i)+B_kalman_model*sum(u_con);
    
    U_con=horzcat(U_con,u_con);
    EPS_VAL=horzcat(EPS_VAL,eps_con);
    
    % Extracting demand for unit time instance for plant simulation
    dc_plant=dc(:,j:j+59);
    j=j+60;
    
    % Non-linear plant simulation with calculated input 'u'
    [qC_all,qT_all,p_bar_all,p_0_all,d_tau]=nonlinear_plant(u_con,dc_plant,Ts,ts,p_0_all,...
        eq2_sym,rhs_eq1_sym,rhs_eq3_sym,dp_sym,dc_sym,d_tau_sym,qc_sym,p_0_sym,tau);
    
    p_0_all_noise=p_0_all+randn(1,length(p_0_all))*std_meas_noise;
    
    
    % Saving variables
    Q_C=horzcat(Q_C,qC_all);
    Q_T=horzcat(Q_T,qT_all);
    P_BAR=horzcat(P_BAR,p_bar_all);
    P_val_0=horzcat(P_val_0,p_0_all(:,2:end));
    D_TAU=horzcat(D_TAU,d_tau);
    P_val_0_noise=horzcat(P_val_0_noise,p_0_all_noise(:,2:end));
    
    % Specifing initial condition for next time instance
    p_0_all=zeros(1,(Ts*60/ts)+1);
    p_0_all(:,1)=P_val_0(:,end);
end

U_con(:,end+1)=U_con(:,end);
H_TANK=P_val_0*1e5/(rho_fluid*g);
H_TANK_noise=P_val_0_noise*1e5/(rho_fluid*g);
H_TANK_esti=P_val_esti*1e5/(rho_fluid*g);
TIME=t(:,1:(T_sec/ts));
TIME_HP=0:T_min;
C_HP=C_con(:,1:T_min+1);
DC_HP=dc(:,1:(T_sec/ts));

D_esti=Demand_C*x_hat(:,1:end-1);
time_plot=linspace(0,T_min-1,length(D_esti));
d2_esti=v1*D_esti;
d5_esti=v2*D_esti;

save('data_control_predictor','U_con','Q_C','Q_T','P_BAR','P_val_0',...
'DC_HP','H_TANK','D_TAU','TIME','TIME_HP','C_HP','T_min','H_TANK_noise',...
'H_TANK_esti','D_esti','time_plot','P_val_0_noise','P_val_esti','d2_esti','d5_esti');
beep

run_time=toc/60;
disp(['Time taken to run MPC control simulation: ',num2str(run_time), ' min']);
disp(' ');

% Plot the results
for_plotting_control_predictor

