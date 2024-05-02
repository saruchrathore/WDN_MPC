clc;
clear all;
close all;

load data_control_predictor


u1_MPC=[];
u2_MPC=[];
C_MPC=[];
for i=1:T_min
u1_MPC=[u1_MPC,U_con(1,i)*ones(1,60)];
u2_MPC=[u2_MPC,U_con(2,i)*ones(1,60)];
C_MPC=[C_MPC,C_HP(i)*ones(1,60)];
end
p1_MPC=P_BAR(1,:);
p2_MPC=P_BAR(5,:);

u1_in=find(u1_MPC<0.02);
u1_MPC(u1_in)=0;

u2_in=find(u2_MPC<0.02);
u2_MPC(u2_in)=0;

eta_p=0.6;                          % Efficiency of pumps
eta_m=0.9;                          % Efficiency of motors
k_eta=1/(eta_p*eta_m*1e3*1e5*3600); % Constant based on efficiency
Ts=1/60;

Cost_pump_1=u1_MPC.*C_MPC.*p1_MPC*k_eta*Ts;
Cost_pump_2=u2_MPC.*C_MPC.*p2_MPC*k_eta*Ts;

j=1;
C_p_1=zeros(1,T_min);
C_p_2=zeros(1,T_min);
for i=1:T_min
    C_p_1(i)=sum(Cost_pump_1(j:j+59));
    C_p_2(i)=sum(Cost_pump_2(j:j+59));
    j=j+60;
end

C_1=zeros(1,T_min);
C_2=zeros(1,T_min);
for i=1:T_min
   C_1(i)= sum(C_p_1(1:i));
   C_2(i)= sum(C_p_2(1:i));
end

C_contr_MPC=[0,C_1+C_2];
C_last_MPC=C_contr_MPC(end)




clearvars -except C_last_MPC C_contr_MPC

load data_control_on_off


u1_MPC=[];
u2_MPC=[];
C_MPC=[];
for i=1:T_min
u1_MPC=[u1_MPC,U_con(1,i)*ones(1,60)];
u2_MPC=[u2_MPC,U_con(2,i)*ones(1,60)];
C_MPC=[C_MPC,C_HP(i)*ones(1,60)];
end
p1_MPC=P_BAR(1,:);
p2_MPC=P_BAR(5,:);

u1_in=find(u1_MPC<0.02);
u1_MPC(u1_in)=0;

u2_in=find(u2_MPC<0.02);
u2_MPC(u2_in)=0;

eta_p=0.6;                          % Efficiency of pumps
eta_m=0.9;                          % Efficiency of motors
k_eta=1/(eta_p*eta_m*1e3*1e5*3600); % Constant based on efficiency
Ts=1/60;

Cost_pump_1=u1_MPC.*C_MPC.*p1_MPC.*k_eta*Ts;
Cost_pump_2=u2_MPC.*C_MPC.*p2_MPC.*k_eta*Ts;

j=1;
C_p_1=zeros(1,T_min);
C_p_2=zeros(1,T_min);
for i=1:T_min
    C_p_1(i)=sum(Cost_pump_1(j:j+59));
    C_p_2(i)=sum(Cost_pump_2(j:j+59));
    j=j+60;
end

C_1=zeros(1,T_min);
C_2=zeros(1,T_min);
for i=1:T_min
   C_1(i)= sum(C_p_1(1:i));
   C_2(i)= sum(C_p_2(1:i));
end

C_contr=[0,C_1+C_2];
C_last_on_off=C_contr(end);

saving_MPC=(C_last_on_off-C_last_MPC)/C_last_on_off*100


set(0, 'DefaultLineLineWidth',2);
set(0, 'DefaultaxesLineWidth',1);
set(0, 'DefaultaxesFontSize',12);
set(0, 'DefaultTextFontSize',12);
set(0, 'DefaultAxesFontName','Times');


t=0:T_min;





figure
plot (t,C_contr_MPC,'g',t,C_contr,'b')
xlabel('Time [min]');ylabel('Accumulated cost of electricity [DKK]')
lgd=legend('NMPC','On/off');
title('Comparison of operational cost between NMPC and On/off control');
xlim([0 T_min]);
% xticks(0:4:T_min)
grid

