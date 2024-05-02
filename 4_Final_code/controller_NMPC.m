%% Controller

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


dc=[d2;d5];  % [m^3/h] Flow demand at consumer node

% Generating demand curve for MPC by calculating average demand in each
% minute
j=1;
dc_con=zeros(2,T_min+Hp);

for i=1:T_min+Hp
    dc_con(:,i)=mean(dc(:,j:j+59),2);
    j=j+60;
end

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

disp(' ');
U_con=[];
Q_C=[];
Q_T=[];
P_BAR=[];
D_TAU=[];
P_val_0=p_0_all(:,1);
EPS_VAL=[];

NMPC_problem


j=1;

for i=1:T_min
    i
    
    % Extracting cost and demand for next Hp instances
    C_Hp=C_con(:,i:i+Hp-1);
    dc_Hp=dc_con(:,i:i+Hp-1);
    
    opti.set_value(p_0_par,p_0_all(:,1));
    opti.set_value(dc_par,dc_Hp);
    opti.set_value(C_par,C_Hp);
    
     
    sol = opti.solve();
    u_all= sol.value(u_var);
    eps_all=sol.value(eps_var);
   
    u_con=u_all(:,1);
    eps_con=eps_all(:,1);
    
    U_con=horzcat(U_con,u_con);
    EPS_VAL=horzcat(EPS_VAL,eps_con);
    
    % Extracting demand for unit time instance for plant simulation
    dc_plant=dc(:,j:j+59);
    j=j+60;
    
    % Non-linear plant simulation with calculated input 'u'
    [qC_all,qT_all,p_bar_all,p_0_all,d_tau]=nonlinear_plant(u_con,dc_plant,Ts,ts,p_0_all,...
        eq2_sym,rhs_eq1_sym,rhs_eq3_sym,dp_sym,dc_sym,d_tau_sym,qc_sym,p_0_sym,tau);
    
    % Saving variables
    Q_C=horzcat(Q_C,qC_all);
    Q_T=horzcat(Q_T,qT_all);
    P_BAR=horzcat(P_BAR,p_bar_all);
    P_val_0=horzcat(P_val_0,p_0_all(:,2:end));
    D_TAU=horzcat(D_TAU,d_tau);
    
    % Specifing initial condition for next time instance
    p_0_all=zeros(1,(Ts*60/ts)+1);
    p_0_all(:,1)=P_val_0(:,end);
end

U_con(:,end+1)=U_con(:,end);
H_TANK=P_val_0*1e5/(rho_fluid*g);
TIME=t(:,1:(T_sec/ts));
TIME_HP=0:T_min;
C_HP=C_con(:,1:T_min+1);
DC_HP=dc(:,1:(T_sec/ts));

save('data_control','U_con','Q_C','Q_T','P_BAR','P_val_0','DC_HP','H_TANK','D_TAU','TIME','TIME_HP','C_HP','T_min');
beep

run_time=toc/60;
disp(['Time taken to run MPC control simulation: ',num2str(run_time), ' min']);
disp(' ');

% Plot the results
for_plotting_control