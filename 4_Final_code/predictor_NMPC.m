y=p_0_noise;


[K,Pc,Pp_kalman]=kfgain_cal(phi_kalman_model,C_kalman_model,Pp_kalman,Q_kalman,R_kalman,n_states_signal+1);        %Calculate Kalman Gain and Covariance matrices

K_kf(:,:,i)=K;                            % Save Kalman Gain for each step
Pc_kf(:,:,i)=Pc;                          % Save Corrected state covariance matrix

%Innovation Process
e(:,i)=y-C_kalman_model*x_hat(:,i);

%Corrected state estimate
x_hat(:,i)=x_hat(:,i)+K*e(:,i);


% Demand_C=[C_signal 0];
% 
% D_esti=Demand_C*x_hat(:,1:end-1);
x_pred=zeros(6,24);
x_pred(:,1)=x_hat(:,i);

for a=1:24
    x_pred(:,a+1)= phi_kalman_model*x_pred(:,a);
end

D_predict=Demand_C*x_pred(:,2:end);
p_esti=x_hat(6,i);