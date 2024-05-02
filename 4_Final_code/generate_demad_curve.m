%% Generating demand curve

    ts=1;           % [s] Sampling time
    T_min=24;           % [min] Simulation time. Provide simulation time in multiple of 24 only.   
    
    % Simulation
    
    T_sec=T_min*60;      % [s] Simulation time in seconds
    N=T_sec/ts;         % [] Simulation instances
    
    % Generating demand curve
    
    t=linspace(0,T_min,(N+1));
    d2=(-0.3*(3*sin((t/(3.83))+5)-4*cos(2*t/(3.83))+8)*0.04)-0.02; % [m^3/s] -ve-> flow out of network
    d5=(-0.18*(3*sin((t/(3.83))+5)-4*cos(2*t/(3.83))+8)*0.04)-0.02; % [m^3/s] -ve-> flow out of network
    
    D_demand=(d2+d5); % Sum of the demand curve