% Consumer demand predictor using Kalman filter



%% Model creation for Kalman filter
[n_states_signal, phi_kalman_model, C_kalman_model, B_kalman_model, tau_s, C_signal]=kalman_model(tau);

return
%% Kalman filter

load('model_data_for_kalman')

%% Dividing the data into 2 parts one to plot for the estimator and other to plot for predictor

T_min=216;

d_tau_2=d_tau(:,(T_min*60)+1:end);
dc_2=dc(:,(T_min*60)+1:end);
dp_2=dp(:,(T_min*60)+1:end);
h_tank_all_2=h_tank_all(:,(T_min*60)+1:end);
p_bar_all_2=p_bar_all(:,(T_min*60)+1:end);
p_ref_all_2=p_ref_all(:,(T_min*60)+1:end);

d_tau=d_tau(:,1:(T_min*60));
dc=dc(:,1:(T_min*60));
dp=dp(:,1:(T_min*60));
h_tank_all=h_tank_all(:,1:(T_min*60));
p_bar_all=p_bar_all(:,1:(T_min*60));
p_ref_all=p_ref_all(:,1:(T_min*60));

v1=mean(dc(1,:)./sum(dc));
v2=mean(dc(2,:)./sum(dc));


%% Kalman filter simulation

T_sim=T_min*60/tau_s;     % Simulation time for signal recreation

std_meas_noise=0.0005;

p_ref_all_noise=p_ref_all+randn(1,length(p_ref_all))*std_meas_noise;
h_tank_all_noise=p_ref_all_noise*1e5/rho_fluid/g;


pump_flow=sum(dp);
y=[];
u=[];

for i=1:tau_s:length(dc)
y=[y p_ref_all_noise(i)]; 
u=[u pump_flow(i)];
end

Pp_kalman=100*diag(1*ones(n_states_signal+1,1));               %Autocovariance for the initial estimate
Q_kalman=(std_process_noise^2)*2*diag(ones(n_states_signal+1,1));%std_state_noise^2;                %Covariance of system noise
Q_kalman(end,end)=Q_kalman(end,end)*(tau*tau_s)^2;
R_kalman=std_meas_noise^2;                    %Covariance of measurement noise

x_hat=zeros(6,1);

for k=1:T_sim
    
    [K,Pc,Pp_kalman]=kfgain_cal(phi_kalman_model,C_kalman_model,Pp_kalman,Q_kalman,R_kalman,n_states_signal+1);        %Calculate Kalman Gain and Covariance matrices
    
    K_kf(:,:,k)=K;                            % Save Kalman Gain for each step
    Pc_kf(:,:,k)=Pc;                          % Save Corrected state covariance matrix
    
    %Innovation Process
    e(:,k)=y(:,k)-C_kalman_model*x_hat(:,k);
    
    %Corrected state estimate
    x_hat(:,k)=x_hat(:,k)+K*e(:,k);
    
    %Predict state for next time step
    x_hat(:,k+1)= phi_kalman_model*x_hat(:,k)+B_kalman_model*u(:,k);
end


Demand_C=[C_signal 0];

D_esti=Demand_C*x_hat(:,1:end-1);

d2_esti=v1*D_esti;
d5_esti=v2*D_esti;

D_demand=sum(dc);
t=linspace(0,T_min,(length(D_demand)+1));
t=t(:,1:T_min*60);
time_plot=linspace(0,T_min-1,length(D_esti));

%% Prediction without measurement signals


x_pred=x_hat(:,k);

for k=1:24  
    x_pred(:,k+1)= phi_kalman_model*x_pred(:,k);   
end

D_predict=Demand_C*x_pred;
t_predict=linspace(T_min-1,T_min+23,length(D_predict));

d2_predict=v1*D_predict;
d5_predict=v2*D_predict;


D_demand_2=sum(dc_2);
t_2=linspace(T_min,240,(length(D_demand_2)+1));
t_2=t_2(:,1:24*60);







save('data_kalman_filter','dp','dc','h_tank_all','D_demand',...
't','time_plot','D_esti','p_ref_all_noise','p_ref_all','x_hat','T_min',...
'D_predict','t_predict','D_demand_2','t_2','d2_predict','d5_predict','dc_2',...
'd2_esti','d5_esti');

for_plotting_kalman

