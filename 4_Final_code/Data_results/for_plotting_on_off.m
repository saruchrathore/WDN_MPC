%% For Plotting Controller

clc;
clear;
close all;


load data_control_on_off


set(0, 'DefaultLineLineWidth',2);
set(0, 'DefaultaxesLineWidth',1);
set(0, 'DefaultaxesFontSize',12);
set(0, 'DefaultTextFontSize',12);
set(0, 'DefaultAxesFontName','Times');


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
h4=H_TANK(1:end-1);


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
plot(t,p1,'g',t,p6,'b');
xlabel('Time [min]');ylabel('Pressure [bar]')
lgd=legend('Pump 1: p_1','Pump 2: p_6');
title('Pressure at the nodes connected to pumps');
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
plot(t,p2,'g',t,p5,'b');
xlabel('Time [min]');ylabel('Pressure [bar]')
lgd=legend('Consumer 1: p_2','Consumer 2: p_5');
title('Pressure at the consumer nodes');
xlim([0 T_min]);
grid


figure
subplot(2,1,1)
stairs(t_Hp,C,'r','LineWidth',2);
xlabel('Time [min]');ylabel('Price [DKK/kW min]')
lgd=legend('Price');
title('Price of electricity');
xlim([0 T_min]);
grid

subplot(2,1,2)
plot(t,h4,'g');
xlabel('Time [min]');ylabel('Level [m]')
lgd=legend('Tank Level: h_4');
title('Level of water in the tank');
xlim([0 T_min]);
grid

