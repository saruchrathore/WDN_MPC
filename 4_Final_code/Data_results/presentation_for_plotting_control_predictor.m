%% For Plotting Controller

clc; clear all; close all;

load data_control_predictor

set(0, 'DefaultLineLineWidth',2);
set(0, 'DefaultaxesLineWidth',1);
set(0, 'DefaultaxesFontSize',12);
set(0, 'DefaultTextFontSize',12);
set(0, 'DefaultAxesFontName','Times');


% rho_fluid=997;  % [kg/m^3] Density of water
% g=9.81;      % [m/s^2] Gravitational constant
% P_val_0_noise=(rho_fluid*g)*H_TANK_noise/1e5;
% P_val_esti=(rho_fluid*g)*H_TANK_esti/1e5;



t=TIME;
t_Hp=TIME_HP;

C=C_HP;

d4=D_TAU;

p1=P_BAR(1,:);
p6=P_BAR(5,:);



d1=U_con(1,:);
d6=U_con(2,:);

p2=P_BAR(2,:);
p5=P_BAR(4,:);


d2=DC_HP(1,:);
d5=DC_HP(2,:);


p4=P_val_0(1:end-1);
p4_noise=P_val_0_noise(1:end-1);
p4_esti=P_val_esti;

h4=H_TANK(1:end-1);
h4_noise=H_TANK_noise(1:end-1);
h4_esti=H_TANK_esti;


p3=P_BAR(3,:);



q1=Q_C(1,:);
q4=Q_C(2,:);



q2=Q_T(1,:);
q3=Q_T(2,:);
q5=Q_T(3,:);
q6=Q_T(4,:);
q7=Q_T(5,:);



figure
subplot(2,1,1)
stairs(t_Hp,d1,'g','LineWidth',2);
hold on
stairs(t_Hp,d6,'b','LineWidth',2);
hold off
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend('Pump 1: u_1','Pump 2: u_2');
title('Flow from the pumps');
xlim([0 T_min]);
grid

subplot(2,1,2)
yyaxis left
plot(t,h4,'g');
xlabel('Time [min]');ylabel('Level [m]')
hold on
yyaxis right
stairs(t_Hp,C,'r','LineWidth',2);
ylim([0.7 1.4]);
ylabel('Price [DKK/kW min]')
lgd=legend('Tank Level: h_4','Price of electricity');
title('Level of water in the tank');
xlim([0 T_min]);
grid





figure
subplot(2,1,1)
plot(t,-d2,'g',t,-d5,'b');
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend('Consumer 1: -d_2','Consumer 2: -d_5');
title('Flow to the consumers');
xlim([0 T_min]);
grid

subplot(2,1,2)
stairs(t_Hp,C,'r','LineWidth',2);
xlabel('Time [min]');ylabel('Price [DKK/kW min]')
lgd=legend('Price');
title('Price of electricity');
xlim([0 T_min]);
grid


figure
subplot(2,1,2)
plot(t,-(d2+d5),'b',time_plot,-D_esti,'g:');
xlabel('Time [min]');ylabel('Flow [m^3/h]')
lgd=legend('Actual consumer demand','Estimated consumer demand from Kalman filter');
title('Comparison of actual consumer demand and estimated consumer demand from Kalman filter predictor');
xlim([0 T_min]);
grid

subplot(2,1,1)
plot(t,p4_noise,'b',t,p4,'r',time_plot,p4_esti,'g:');
xlabel('Time [min]');ylabel('Pressure [bar]')
lgd=legend('Measurement with noise','Measurement','Kalman filter estimate');
title('Comparison of actual tank pressure measurement and estimated measurement from Kalman filter predictor');
xlim([0 T_min]);
grid