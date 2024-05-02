%% MPC plant simulation
% Running non-linear plant simulation with input provided by MPC

function [qC_all,qT_all,p_bar_all,p_0_all,d_tau]=nonlinear_plant(u_con,dc_plant,Ts,ts,p_0_all,...
    eq2_sym,rhs_eq1_sym,rhs_eq3_sym,dp_sym,dc_sym,d_tau_sym,qc_sym,p_0_sym,tau)

N=Ts*60/ts;

qC_all=zeros(2,N);
qT_all=zeros(5,N);
p_bar_all=zeros(5,N);
d_tau=-sum(u_con+dc_plant,1);

for i=1:N
    
    eq2_sym_sub=subs(eq2_sym,[dp_sym;dc_sym;d_tau_sym],[u_con;dc_plant(:,i);d_tau(:,i)]);
    
    [q1,q4]=vpasolve(eq2_sym_sub);
    
    qC=[q1;q4];
    
    qC_all(:,i)=qC;
    
    qT=double(subs(rhs_eq1_sym,[qc_sym;dp_sym;dc_sym],[qC;u_con;dc_plant(:,i)]));
    
    qT_all(:,i)=qT;
    
    
    p_bar_all(:,i)=double(subs(rhs_eq3_sym,[qc_sym;dp_sym;dc_sym;p_0_sym],[qC;u_con;dc_plant(:,i);p_0_all(:,i)]));
    
    p_0_all(:,i+1)=p_0_all(:,i)-tau*d_tau(:,i)*ts;
    
end
