%% For Plotting

load('data_kalman_filter')

set(0, 'DefaultLineLineWidth',2);
set(0, 'DefaultaxesLineWidth',1);
set(0, 'DefaultaxesFontSize',12);
set(0, 'DefaultTextFontSize',12);
set(0, 'DefaultAxesFontName','Times');




%% Plot nonlinear simulation

d1=dp(1,:);
d6=dp(2,:);

d2=dc(1,:);
d5=dc(2,:);

d2_2=dc_2(1,:);
d5_2=dc_2(2,:);

p4=p_ref_all;
h4=h_tank_all;


figure
subplot(3,1,1)
plot(t,d1,'g',t,d6,'b');
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend('Pump 1: u_1','Pump 2: u_2');
title('Flow from the pumps');
xlim([0 T_min]);
grid

subplot(3,1,2)
plot(t,-d2,'g',t,-d5,'b');
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend('Consumer 1: -d_2','Consumer 2: -d_5');
title('Flow to the consumers');
xlim([0 T_min]);
grid


subplot(3,1,3)
plot(t,h4,'g--');
xlabel('Time [min]');ylabel('Level [m]')
lgd=legend('Tank Level: h_4');
title('Level of water in the tank');
xlim([0 T_min]);
grid

%% Plot Kalman filter results


p_tank_noise=p_ref_all_noise;
p_tank=p_ref_all;
p_esti=x_hat(6,1:end-1);

figure
plot(t,-D_demand,'b',time_plot,-D_esti,'g:')
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend('Actual consumer demand','Estimated consumer demand from Kalman filter');
title('Comparison of actual consumer demand and estimated consumer demand from Kalman filter');
xlim([0 T_min]);
grid


figure
plot(t,p_tank_noise,'b',t,p_tank,'r',time_plot,p_esti,'g:')
xlabel('Time [min]');ylabel('Pressure [bar]')
lgd=legend('Measurement with noise','Measurement','Kalman filter estimate');
title('Comparison of actual tank pressure measurement and estimated measurement from Kalman filter');
xlim([0 T_min]);
grid



figure
plot_1=plot(t,-D_demand,'b',time_plot,-D_esti,'g:');
hold on
plot_2=plot(t_2,-D_demand_2,'b',t_predict,-D_predict,'m:');
hold off
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend([plot_1(1) plot_1(2) plot_2(2)],'Actual consumer demand','Estimated consumer demand','Predicted consumer demand');
title('Comparison of actual consumer demand and predicted consumer demand without measurements');
xlim([T_min-24*3 T_min+24]);
grid



figure
subplot(2,1,1)
plot_1=plot(t,-d2,'b',time_plot,-d2_esti,'g:');
hold on
plot_2=plot(t_2,-d2_2,'b',t_predict,-d2_predict,'m:');
hold off
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend([plot_1(1) plot_1(2) plot_2(2)],'Actual d2 demand','Estimated d2 demand','Predicted d2 demand');
title('Comparison of d2 consumer demand: actual and predicted  without measurements');
xlim([T_min-24*3 T_min+24]);
grid

subplot(2,1,2)
plot_1=plot(t,-d5,'b',time_plot,-d5_esti,'g:');
hold on
plot_2=plot(t_2,-d5_2,'b',t_predict,-d5_predict,'m:');
hold off
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend([plot_1(1) plot_1(2) plot_2(2)],'Actual d5 demand','Estimated d5 demand','Predicted d5 demand');
title('Comparison of d5 consumer demand: actual and predicted  without measurements');
xlim([T_min-24*3 T_min+24]);
grid



